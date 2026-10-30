"""Verify topology copying and ion regrouping preserve the fitted system."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from vaasu.topology import prepare_topology


def gro_record(index, residue):
    return (
        f"{index:5d}{residue:<5}{residue:>5}{index:5d}"
        f"{index / 10:8.3f}{0.2:8.3f}{0.3:8.3f}\n"
    )


class PrepareTopologyTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.vanilla = self.root / "vanilla"
        self.boresch = self.root / "boresch"
        self.complex = self.root / "complex"
        for directory in (self.vanilla, self.boresch, self.complex):
            directory.mkdir()
        (self.vanilla / "toppar").mkdir()
        (self.vanilla / "toppar" / "ligand.itp").write_text(
            '; nested include\n#include "posre.itp"\n'
        )
        (self.vanilla / "toppar" / "posre.itp").write_text("; position restraints\n")
        force_field = self.vanilla / "example.ff"
        force_field.mkdir()
        (force_field / "forcefield.itp").write_text('#include "ffnonbonded.itp"\n')
        (force_field / "ffnonbonded.itp").write_text("; nonbonded parameters\n")
        self.prefix = (
            '; topology\n#include "./example.ff/forcefield.itp"\n'
            '#include "toppar/ligand.itp"\n[ system ]\nTest system\n'
            "[ molecules ] ; molecule list\n; molecule counts\nProtein 1\nunk 1\nSOL 1\n"
        )
        self.topology = (
            self.prefix
            + "CL 1 ; first chloride\nNA 1\n; ion comment\nCL 1 ; last chloride\n"
        )
        (self.vanilla / "topol.top").write_text(self.topology)
        self.records = [
            gro_record(index, residue)
            for index, residue in enumerate(("GLY", "unk", "SOL", "CL", "NA", "CL"), 1)
        ]
        self.coordinates = "Fitted frame\n6\n" + "".join(self.records) + "2.0 2.0 2.0\n"
        (self.boresch / "ClosestRestraintFrame.gro").write_text(self.coordinates)
        self.restraints = (
            "\n; restraints\n[ intermolecular_interactions ]\n[ bonds ]\n"
            "2 1 6 0.5 0 0.5 100\n"
        )
        (self.boresch / "BoreschRestraint.top").write_text(self.restraints)

    def prepare(self):
        return prepare_topology(self.vanilla, self.boresch, self.complex)

    def test_preserves_anchors_records_box_restraints_and_sources(self):
        source_files = {
            path: path.read_bytes()
            for directory in (self.vanilla, self.boresch)
            for path in directory.rglob("*")
            if path.is_file()
        }
        working_directory = Path.cwd()
        self.assertEqual(self.prepare(), self.complex)
        self.assertEqual(Path.cwd(), working_directory)
        grouped = (self.complex / "complex.gro").read_text().splitlines(keepends=True)
        self.assertEqual(grouped[:5], self.coordinates.splitlines(keepends=True)[:5])
        self.assertEqual(
            grouped[5:-1], [self.records[4], self.records[3], self.records[5]]
        )
        self.assertEqual(grouped[-1], "2.0 2.0 2.0\n")
        self.assertEqual(
            (self.complex / "complex_ungrouped.gro").read_text(), self.coordinates
        )
        self.assertEqual(
            (self.complex / "complex_ungrouped.top").read_text(),
            self.topology + "\n" + self.restraints,
        )
        self.assertEqual(
            (self.complex / "BoreschRestraint.top").read_text(), self.restraints
        )
        topology = (self.complex / "complex.top").read_text()
        self.assertTrue(topology.startswith(self.prefix))
        for comment in ("; first chloride", "; ion comment", "; last chloride"):
            self.assertIn(comment, topology)
        self.assertIn("NA           1\nCL           2\n", topology)
        self.assertTrue(topology.endswith(self.restraints))
        self.assertEqual(topology.count("[ intermolecular_interactions ]"), 1)
        for path, contents in source_files.items():
            self.assertEqual(path.read_bytes(), contents)
        for name in (
            "toppar/ligand.itp",
            "toppar/posre.itp",
            "example.ff/forcefield.itp",
            "example.ff/ffnonbonded.itp",
        ):
            self.assertEqual(
                (self.complex / name).read_bytes(), (self.vanilla / name).read_bytes()
            )

    def test_handles_absent_ion_species_and_no_ions(self):
        for residues in (("GLY", "NA"), ("GLY", "CL"), ("GLY",)):
            with (
                self.subTest(residues=residues),
                TemporaryDirectory(dir=self.root) as destination,
            ):
                records = [
                    gro_record(index, residue)
                    for index, residue in enumerate(residues, 1)
                ]
                coordinates = f"Frame\n{len(records)}\n" + "".join(records) + "2 2 2\n"
                (self.boresch / "ClosestRestraintFrame.gro").write_text(coordinates)
                topology = "[ molecules ]\nProtein 1\n" + "".join(
                    f"{name} 1\n" for name in residues[1:]
                )
                (self.vanilla / "topol.top").write_text(topology)
                prepare_topology(self.vanilla, self.boresch, destination)
                self.assertEqual(
                    (Path(destination) / "complex.gro").read_text(), coordinates
                )

    def test_copies_other_direct_local_includes(self):
        (self.vanilla / "parameters").mkdir()
        (self.vanilla / "parameters" / "extra.itp").write_text("; extra parameters\n")
        (self.vanilla / "topol.top").write_text(
            "#include <parameters/extra.itp>\n" + self.topology
        )
        self.prepare()
        self.assertEqual(
            (self.complex / "parameters" / "extra.itp").read_text(),
            "; extra parameters\n",
        )

    def test_rejects_mismatched_counts_before_copying(self):
        (self.vanilla / "topol.top").write_text(self.topology.replace("NA 1", "NA 2"))
        with self.assertRaisesRegex(ValueError, "Ion counts disagree"):
            self.prepare()
        self.assertEqual(list(self.complex.iterdir()), [])

    def test_rejects_ions_before_non_ion_atoms(self):
        records = self.records.copy()
        records[0], records[3] = records[3], records[0]
        (self.boresch / "ClosestRestraintFrame.gro").write_text(
            "Frame\n6\n" + "".join(records) + "2 2 2\n"
        )
        with self.assertRaisesRegex(ValueError, "preserve anchor indices"):
            self.prepare()
        self.assertEqual(list(self.complex.iterdir()), [])

    def test_rejects_ions_before_non_ion_molecules(self):
        (self.vanilla / "topol.top").write_text(
            "[ molecules ]\nCL 2\nProtein 1\nNA 1\n"
        )
        with self.assertRaisesRegex(ValueError, "non-ion molecules"):
            self.prepare()
        self.assertEqual(list(self.complex.iterdir()), [])

    def test_preserves_external_include_targets_and_sources(self):
        external = self.root / "extra.itp"
        external.write_text("; external parameters\n")
        (self.vanilla / "topol.top").write_text(
            f'#include "../extra.itp"\n#include "{external}"\n' + self.topology
        )
        self.prepare()
        topology = (self.complex / "complex.top").read_text()
        self.assertEqual(topology.count(f'#include "{external}"'), 2)
        self.assertEqual(external.read_text(), "; external parameters\n")
        self.assertFalse((self.complex / "extra.itp").exists())

    def test_preserves_existing_restraints_after_grouped_ions(self):
        (self.vanilla / "topol.top").write_text(self.topology + self.restraints)
        self.prepare()
        topology = (self.complex / "complex.top").read_text()
        self.assertLess(
            topology.index("NA           1\nCL           2\n"),
            topology.index(self.restraints),
        )
        self.assertTrue(topology.endswith(self.restraints))
        self.assertEqual(topology.count(self.restraints), 2)

    def test_rejects_missing_or_duplicate_molecule_sections(self):
        for topology in (
            "[ system ]\nSystem\n",
            self.topology + "[ molecules ]\nNA 1\n",
        ):
            with self.subTest(topology=topology):
                (self.vanilla / "topol.top").write_text(topology)
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assertEqual(list(self.complex.iterdir()), [])

    def test_rejects_malformed_atom_records(self):
        for coordinates in (
            self.coordinates.replace(self.records[0], "bad atom record\n"),
            self.coordinates.replace(self.records[0][20:28], "     nan"),
        ):
            with self.subTest(coordinates=coordinates):
                (self.boresch / "ClosestRestraintFrame.gro").write_text(coordinates)
                with self.assertRaisesRegex(ValueError, "Malformed GRO"):
                    self.prepare()
                self.assertEqual(list(self.complex.iterdir()), [])

    def test_rejects_missing_restraint_outputs(self):
        for name in ("ClosestRestraintFrame.gro", "BoreschRestraint.top"):
            with self.subTest(name=name):
                path = self.boresch / name
                contents = path.read_bytes()
                path.unlink()
                with self.assertRaisesRegex(FileNotFoundError, name):
                    self.prepare()
                self.assertEqual(list(self.complex.iterdir()), [])
                path.write_bytes(contents)

    def test_rejects_existing_files(self):
        path = self.complex / "complex.top"
        path.write_text("existing data")
        with self.assertRaises(FileExistsError):
            self.prepare()
        self.assertEqual(path.read_text(), "existing data")


if __name__ == "__main__":
    unittest.main()
