"""Verify alchemical parameters, counterion placement, and output safeguards."""

import shutil
import unittest
from contextlib import redirect_stdout
from io import StringIO
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import MDAnalysis as mda

from vaasu.alchemy import prepare_alchemical_system
from vaasu.windows import create_windows


def gro_record(number, residue_number, residue, atom, x, *, precision=3):
    width = precision + 5
    return (
        f"{residue_number:5d}{residue:<5}{atom:>5}{number:5d}"
        f"{x:{width}.{precision}f}{1.0:{width}.{precision}f}{1.0:{width}.{precision}f}\n"
    )


def section_rows(text, name):
    """Read numeric output for assertions independently of the implementation."""
    rows = []
    section = ""
    for line in text.splitlines():
        content = line.split(";")[0].strip()
        if content.startswith("["):
            section = content.strip("[] ")
        elif content and section == name:
            rows.append(content.split())
    return rows


class AlchemyTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.toppar = self.root / "toppar"
        self.toppar.mkdir()
        self.ligand = self.toppar / "ligand.itp"
        self.atomtypes = self.toppar / "ligand_atomtypes.itp"
        self.bonded = "\n[ bonds ]\n; keep this bonded term\n1 2 1 0.15 1000\n"
        self.write_ligand(1)
        self.atomtypes.write_text(
            "[ atomtypes ]\n; name number mass charge ptype sigma epsilon\n"
            "c1 6 12.0 0 A 0.3 0.4\n"
            "n1 7 14.0 0 A 0.4 0.9 ; nitrogen\n"
        )
        (self.root / "forcefield.itp").write_text(
            "[ defaults ]\n1 2 yes 0.5 0.8333\n"
            '#include "nested.itp"\n'
            "[ atomtypes ]\nNa 11 22.99 0 A 0.3 0.1\nCl 17 35.45 0 A 0.4 0.4\n"
        )
        self.restraints = (
            "; fitted Boresch restraints\n[ intermolecular_interactions ]\n"
            "[ bonds ]\n2 1 6 0.5 0 0.5 100\n"
            "[ angles ]\n1 2 3 1 90 0 90 20\n"
        )
        self.topology = (
            '#include "forcefield.itp"\n'
            '#include "./toppar/ligand_atomtypes.itp"\n'
            '#include "./toppar/ligand.itp"\n'
            "; toppar/ligand.itp in a comment must remain unchanged\n"
            "[ system ]\nFixture\n[ molecules ]\nProtein 1\nunk 1\nSOL 1\n"
            "NA 3\nCL 3 ; ordinary chlorides\n\n" + self.restraints
        )
        (self.root / "complex.top").write_text(self.topology)
        specifications = [
            (1, "GLY", "CA", 0.1),
            (2, "unk", "C1", 1.0),
            (2, "unk", "N1", 1.1),
            (3, "SOL", "OW", 0.5),
            (4, "NA", "NA", 1.7),
            (5, "NA", "NA", 1.8),
            (6, "NA", "NA", 1.9),
            (7, "CL", "CL", 1.7),
            (8, "CL", "CL", 1.8),
            (9, "CL", "CL", 1.9),
        ]
        self.records = [
            gro_record(index, *specification)
            for index, specification in enumerate(specifications, 1)
        ]
        self.coordinates = "Fixture\n10\n" + "".join(self.records) + "2.0 2.0 2.0\n"
        (self.root / "complex.gro").write_text(self.coordinates)

    def write_ligand(self, charge):
        self.ligand.write_text(
            "[ moleculetype ]\nunk 3\n\n[ atoms ]\n"
            f"1 c1 1 unk C1 1 {charge - 0.25} 12.0 ; carbon\n"
            "2 n1 1 unk N1 1 0.25 14.0\n" + self.bonded
        )

    def prepare(self, **kwargs):
        return prepare_alchemical_system(
            self.root,
            ligand_selection=kwargs.get("ligand_selection", "resname unk"),
            ligand_itp=kwargs.get("ligand_itp", "toppar/ligand.itp"),
            atomtype_itp=kwargs.get("atomtype_itp", "toppar/ligand_atomtypes.itp"),
        )

    def assert_no_outputs(self):
        self.assertFalse(list(self.root.glob("complex_*.gro")))
        self.assertFalse(list(self.root.glob("complex_*.top")))
        self.assertFalse(list(self.toppar.glob("*_alchemical.itp")))

    def test_writes_parameters_and_preserves_sources_restraints_and_cwd(self):
        sources = {
            path: path.read_bytes() for path in self.root.rglob("*") if path.is_file()
        }
        working_directory = Path.cwd()
        self.assertEqual(self.prepare(), self.root)
        self.assertEqual(Path.cwd(), working_directory)
        for path, content in sources.items():
            self.assertEqual(path.read_bytes(), content)
        ligand_text = (self.toppar / "ligand_alchemical.itp").read_text()
        self.assertTrue(ligand_text.endswith(self.bonded))
        atoms = section_rows(ligand_text, "atoms")
        self.assertEqual(
            atoms[0],
            ["1", "c1", "1", "unk", "C1", "1", "0.75", "12.0", "dc1", "0.0", "12.0"],
        )
        self.assertEqual(atoms[1][-3:], ["dn1", "0.0", "14.0"])
        atomtypes_text = (self.toppar / "ligand_atomtypes_alchemical.itp").read_text()
        types = {row[0]: row for row in section_rows(atomtypes_text, "atomtypes")}
        self.assertEqual(types["dc1"], ["dc1", "6", "12.0", "0", "A", "0.0", "0.0"])
        pairs = section_rows(atomtypes_text, "nonbond_params")
        self.assertEqual(len(pairs), 6)
        cross = next(row for row in pairs if row[:2] == ["dc1", "dn1"])
        self.assertAlmostEqual(float(cross[3]), 0.35)
        self.assertAlmostEqual(float(cross[4]), 0.6)
        topology = (self.root / "complex_alchemical.top").read_text()
        self.assertTrue(topology.endswith(self.restraints))
        self.assertEqual(topology.count("[ intermolecular_interactions ]"), 1)
        self.assertIn('#include "toppar/ligand_alchemical.itp"', topology)
        self.assertIn('#include "toppar/ligand_atomtypes_alchemical.itp"', topology)
        self.assertIn(
            "; toppar/ligand.itp in a comment must remain unchanged", topology
        )
        self.assertIn("CL 2 ; ordinary chlorides", topology)
        self.assertIn("1 1 1000 1000 1000", topology)
        self.assertEqual(
            {path.name for path in self.toppar.glob("*_alchemical.itp")},
            {"ligand_alchemical.itp", "ligand_atomtypes_alchemical.itp"},
        )

    def test_neutral_and_signed_charges_preserve_counterion_placement(self):
        for charge in (0, 1, -1, 2, -2, 3, -3):
            with self.subTest(charge=charge), TemporaryDirectory() as destination:
                self.write_ligand(charge)
                root = Path(destination) / "system"
                shutil.copytree(self.root, root)
                prepare_alchemical_system(
                    root,
                    ligand_selection="resname unk",
                    ligand_itp="toppar/ligand.itp",
                    atomtype_itp="toppar/ligand_atomtypes.itp",
                )
                path = root / "complex_alchemical.gro"
                topology = (root / "complex_alchemical.top").read_text()
                self.assertTrue(topology.endswith(self.restraints))
                if charge == 0:
                    self.assertEqual(path.read_text(), self.coordinates)
                    self.assertNotIn("B_CL", topology)
                    self.assertNotIn("B_NA", topology)
                    continue
                name, ordinary, ion_charge, mass = (
                    ("B_CL", "CL", -1, 35.45)
                    if charge > 0
                    else ("B_NA", "NA", 1, 22.99)
                )
                universe = mda.Universe(str(path), to_guess=())
                try:
                    alchemical = universe.select_atoms(f"resname {name}")
                    self.assertEqual(alchemical.n_atoms, abs(charge))
                    self.assertEqual(
                        alchemical.indices.tolist(), list(range(10 - abs(charge), 10))
                    )
                    actual_positions = sorted(
                        tuple(position) for position in alchemical.positions
                    )
                    expected_positions = sorted(
                        [(0, 0, 0), (10, 0, 0), (0, 10, 0)][: abs(charge)]
                    )
                    self.assertEqual(actual_positions, expected_positions)
                    self.assertEqual(
                        universe.select_atoms(f"resname {ordinary}").n_atoms,
                        3 - abs(charge),
                    )
                finally:
                    universe.trajectory.close()
                molecules = dict(section_rows(topology, "molecules"))
                self.assertEqual(int(molecules[name]), abs(charge))
                self.assertEqual(int(molecules[ordinary]), 3 - abs(charge))
                ion_atom = section_rows(topology, "atoms")[0]
                self.assertEqual(float(ion_atom[6]), ion_charge)
                self.assertEqual(float(ion_atom[9]), 0)
                self.assertEqual(float(ion_atom[7]), mass)
                self.assertEqual(ion_atom[1], ion_atom[8])
                self.assertEqual(
                    path.read_text().splitlines(keepends=True)[2:6], self.records[:4]
                )

    def test_charge_is_conserved_at_every_global_lambda_state(self):
        for charge in (0, 1, -1, 2, -2, 3, -3):
            with self.subTest(charge=charge), TemporaryDirectory() as destination:
                self.write_ligand(charge)
                root = Path(destination) / "system"
                shutil.copytree(self.root, root)
                prepare_alchemical_system(
                    root,
                    ligand_selection="resname unk",
                    ligand_itp="toppar/ligand.itp",
                    atomtype_itp="toppar/ligand_atomtypes.itp",
                )
                ligand_atoms = section_rows(
                    (root / "toppar/ligand_alchemical.itp").read_text(), "atoms"
                )
                ligand_a = fsum(float(atom[6]) for atom in ligand_atoms)
                ligand_b = fsum(float(atom[9]) for atom in ligand_atoms)
                ion_atoms = section_rows(
                    (root / "complex_alchemical.top").read_text(), "atoms"
                )
                ions_a = abs(charge) * float(ion_atoms[0][6]) if charge else 0
                ions_b = abs(charge) * float(ion_atoms[0][9]) if charge else 0
                self.assertAlmostEqual(ligand_a + ions_a, 0)
                self.assertAlmostEqual(ligand_b + ions_b, 0)
                with redirect_stdout(StringIO()):
                    windows = create_windows(root, temperature=298.15)
                for state_index, window in enumerate(windows):
                    mdp = next((window / "prod").glob("*.mdp")).read_text()
                    coul_line = next(
                        line
                        for line in mdp.splitlines()
                        if line.startswith("coul-lambdas")
                    )
                    coul = float(coul_line.split("=", 1)[1].split()[state_index])
                    ligand_charge = (1 - coul) * ligand_a + coul * ligand_b
                    ion_charge = (1 - coul) * ions_a + coul * ions_b
                    self.assertAlmostEqual(ligand_charge + ion_charge, 0)
                    if window.name.startswith("rest."):
                        self.assertAlmostEqual(ligand_charge, charge)
                    elif window.name.startswith("vdw."):
                        self.assertAlmostEqual(ligand_charge, 0)

    def test_selects_farthest_ions_and_breaks_ties_by_original_order(self):
        records = self.records.copy()
        records[7] = gro_record(8, 7, "CL", "CL", 1.9)
        (self.root / "complex.gro").write_text(
            "Fixture\n10\n" + "".join(records) + "2 2 2\n"
        )
        self.prepare()
        output = (self.root / "complex_alchemical.gro").read_text().splitlines()
        # Both chlorides 8 and 10 are equally far away; original atom 8 wins.
        self.assertEqual(int(output[-2][15:20]), 8)
        self.assertEqual(output[-2][5:10].strip(), "B_CL")
        self.assertEqual(output[2:6], [line.rstrip("\n") for line in self.records[:4]])

    def test_preserves_precision_velocities_box_and_nonselected_records(self):
        records = []
        for record in self.records:
            records.append(
                gro_record(
                    int(record[15:20]),
                    int(record[:5]),
                    record[5:10].strip(),
                    record[10:15].strip(),
                    float(record[20:28]),
                    precision=5,
                ).rstrip("\n")
                + "  0.123456  0.234567  0.345678\n"
            )
        box = "2.0 2.0 2.0 0 0 0.2 0 0.3 0.4\n"
        (self.root / "complex.gro").write_text("Fixture\n10\n" + "".join(records) + box)
        self.prepare()
        output = (
            (self.root / "complex_alchemical.gro").read_text().splitlines(keepends=True)
        )
        self.assertEqual(output[2:-2], records[:-1])
        self.assertEqual(output[-1], box)
        self.assertEqual(output[-2][50:], records[-1][50:])
        self.assertEqual(output[-2][20:50], f"{0:10.5f}{0:10.5f}{0:10.5f}")

    def test_rejects_missing_inputs_without_outputs(self):
        for path in (
            self.ligand,
            self.atomtypes,
            self.root / "complex.top",
            self.root / "complex.gro",
        ):
            with self.subTest(path=path):
                original = path.read_bytes()
                path.unlink()
                with self.assertRaises(FileNotFoundError):
                    self.prepare()
                self.assert_no_outputs()
                path.write_bytes(original)

    def test_rejects_fractional_charge_or_missing_counterions_without_outputs(self):
        for charge, message in ((0.5, "integer"), (4, "counterions")):
            with self.subTest(charge=charge):
                self.write_ligand(charge)
                with self.assertRaisesRegex(ValueError, message):
                    self.prepare()
                self.assert_no_outputs()

    def test_rejects_output_conflicts_including_dangling_symlinks(self):
        path = self.toppar / "ligand_alchemical.itp"
        path.write_text("keep previous result")
        with self.assertRaises(FileExistsError):
            self.prepare()
        self.assertEqual(path.read_text(), "keep previous result")
        self.assertFalse((self.root / "complex_alchemical.gro").exists())
        path.unlink()
        path.symlink_to("missing.itp")
        with self.assertRaises(FileExistsError):
            self.prepare()
        self.assertTrue(path.is_symlink())
        self.assertFalse((self.root / "complex_alchemical.gro").exists())

    def test_rejects_missing_or_duplicate_restraints(self):
        variants = [
            self.topology.replace(self.restraints, ""),
            self.topology + "\n[ intermolecular_interactions ]\n",
        ]
        for text in variants:
            with self.subTest(text=text):
                (self.root / "complex.top").write_text(text)
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assert_no_outputs()

    def test_neutral_ligands_still_require_valid_system_sections(self):
        self.write_ligand(0)
        (self.root / "complex.top").write_text(
            self.topology.replace(self.restraints, "")
        )
        with self.assertRaisesRegex(ValueError, "intermolecular_interactions"):
            self.prepare()
        self.assert_no_outputs()

    def test_io_failure_preserves_sources_and_partial_outputs(self):
        original_open = Path.open

        def fail_later(path, *args, **kwargs):
            if path.name == "ligand_alchemical.itp":
                raise OSError("disk write failed")
            return original_open(path, *args, **kwargs)

        source = (self.root / "complex.top").read_bytes()
        with (
            patch.object(Path, "open", fail_later),
            self.assertRaisesRegex(OSError, "disk write failed"),
        ):
            self.prepare()
        self.assertTrue((self.root / "complex_alchemical.top").is_file())
        self.assertEqual((self.root / "complex.top").read_bytes(), source)


if __name__ == "__main__":
    unittest.main()
