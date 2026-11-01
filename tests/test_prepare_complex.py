"""Check setup safeguards and the restraint-generation integration."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import Mock, call, patch

from vaasu.boresch import generate_boresch_restraints
from vaasu.complex import prepare_complex


class PrepareComplexTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.vanilla = self.root / "vanilla"
        self.vanilla.mkdir()
        for name in ("prod.tpr", "prod.xtc"):
            (self.vanilla / name).touch()
        (self.vanilla / "toppar").mkdir()
        (self.vanilla / "toppar" / "ligand.itp").write_text(
            "[ moleculetype ]\nunk 3\n[ atoms ]\n1 c1 1 unk unk 1 0 12.0\n"
        )
        (self.vanilla / "toppar" / "ligand_atomtypes.itp").write_text(
            "[ atomtypes ]\nc1 6 12.0 0 A 0.3 0.4\n"
        )
        force_field = self.vanilla / "example.ff"
        force_field.mkdir()
        (force_field / "forcefield.itp").write_text(
            "[ defaults ]\n1 2 yes 0.5 0.8333\n"
        )
        (self.vanilla / "topol.top").write_text(
            '#include "./example.ff/forcefield.itp"\n'
            '#include "toppar/ligand_atomtypes.itp"\n'
            '#include "toppar/ligand.itp"\n'
            "[ system ]\nFixture\n"
            "[ molecules ]\nProtein 1\nunk 1\nCL 1\nNA 1\nCL 1\n"
        )
        self.abfe = self.root / "nested" / "abfe"
        self.generator = self.enterContext(
            patch("vaasu.complex.generate_boresch_restraints")
        )
        self.generator.side_effect = self.write_restraints

    def write_restraints(self, topology, trajectory, output_dir, **kwargs):
        output_dir.mkdir()
        records = [
            f"{index:5d}{name:<5}{name:>5}{index:5d}{0.1:8.3f}{0.2:8.3f}{0.3:8.3f}\n"
            for index, name in enumerate(("GLY", "unk", "CL", "NA", "CL"), 1)
        ]
        (output_dir / "ClosestRestraintFrame.gro").write_text(
            "Fitted frame\n5\n" + "".join(records) + "2.0 2.0 2.0\n"
        )
        (output_dir / "BoreschRestraint.top").write_text(
            "; restraints\n[ intermolecular_interactions ]\n[ bonds ]\n1 2 6 0.5 0 0.5 100\n"
        )

    def prepare(self, **kwargs):
        return prepare_complex(
            self.vanilla,
            self.abfe,
            ligand_selection="resname unk",
            temperature=kwargs.get("temperature", 298.15),
            ligand_itp="toppar/ligand.itp",
            atomtype_itp="toppar/ligand_atomtypes.itp",
            **{key: value for key, value in kwargs.items() if key != "temperature"},
        )

    def test_prepares_complex_and_passes_explicit_paths(self):
        working_directory = Path.cwd()
        result = self.prepare()

        self.assertEqual(result, self.abfe / "complex")
        self.assertEqual(
            {path.name for path in result.iterdir()},
            {
                "job_complex_archer_test.sh",
                "job_complex_archer.sh",
                "job_complex_archer_contd.sh",
                "job_complex_archer_lowpriority.sh",
                "job_complex_archer_contd_lowpriority.sh",
                "plumed.dat",
                "index_alchemical.ndx",
                "simulations_list.txt",
                *[f"rest.{number:02d}" for number in range(12)],
                *[f"coul.{number:02d}" for number in range(11)],
                *[f"vdw.{number:02d}" for number in range(21)],
                "complex_alchemical.gro",
                "complex_alchemical.top",
                "complex.gro",
                "complex.top",
                "complex_ungrouped.gro",
                "complex_ungrouped.top",
                "BoreschRestraint.top",
                "toppar",
                "example.ff",
            },
        )
        self.assertIn(
            "NA           1\nCL           2\n", (result / "complex.top").read_text()
        )
        self.assertTrue((result / "toppar" / "ligand.itp").is_file())
        self.assertTrue((result / "example.ff" / "forcefield.itp").is_file())
        self.assertEqual(Path.cwd(), working_directory)
        self.generator.assert_called_once_with(
            self.vanilla / "prod.tpr",
            self.vanilla / "prod.xtc",
            self.abfe / "boresch",
            ligand_selection="resname unk",
            temperature=298.15,
        )

    def test_forwards_custom_coupling_and_window_settings(self):
        result = self.prepare(
            temperature=310.0,
            job_name="DOR_10",
            solvent_selection="resname CL",
            lambda_schedule={"rest": [0.0, 1.0], "coul": [0.0, 1.0], "vdw": [0.0, 1.0]},
        )
        self.assertEqual(
            len((result / "simulations_list.txt").read_text().splitlines()), 6
        )
        job = (result / "job_complex_archer.sh").read_text()
        self.assertIn("#SBATCH --job-name=DOR_10\n", job)
        self.assertIn("#SBATCH --ntasks=6\n", job)
        self.assertIn("#SBATCH --nodes=3\n", job)
        self.assertIn("windows=(rest.00 rest.01 coul.00 coul.01 vdw.00 vdw.01)", job)

        index = (result / "index_alchemical.ndx").read_text()
        self.assertIn("[ Solvent ]\n4 5", index)
        self.assertIn("[ Solute ]\n1 2 3", index)
        mdp = (result / "coul.00/prod/coul.prod.00.mdp").read_text()
        self.assertIn("init-lambda-state        = 2", mdp)
        temperature_line = next(
            line for line in mdp.splitlines() if line.startswith("ref-t")
        )
        self.assertEqual(temperature_line.split("=", 1)[1].split(), ["310.0", "310.0"])

    def test_rejects_existing_outputs_without_touching_them(self):
        for name in ("complex", "boresch"):
            with self.subTest(name=name):
                destination = self.abfe / name
                destination.mkdir(parents=True)
                sentinel = destination / "keep.txt"
                sentinel.write_text("existing data")
                with self.assertRaises(FileExistsError):
                    self.prepare()
                self.assertEqual(sentinel.read_text(), "existing data")
                self.generator.assert_not_called()
                sentinel.unlink()
                destination.rmdir()

    def test_missing_inputs_do_not_create_abfe_directory(self):
        for name in (
            "prod.tpr",
            "prod.xtc",
            "topol.top",
            "example.ff/forcefield.itp",
            "toppar/ligand.itp",
            "toppar/ligand_atomtypes.itp",
        ):
            with self.subTest(name=name):
                missing = self.vanilla / name
                contents = missing.read_bytes()
                missing.unlink()
                with self.assertRaisesRegex(FileNotFoundError, name):
                    self.prepare()
                self.assertFalse(self.abfe.exists())
                self.generator.assert_not_called()
                missing.write_bytes(contents)

    def test_alchemical_parameters_are_not_parsed_before_fitting(self):
        self.generator.side_effect = RuntimeError("fitting reached")
        with (
            patch("vaasu.alchemy._read_parameters") as read_parameters,
            self.assertRaisesRegex(RuntimeError, "fitting reached"),
        ):
            self.prepare()
        self.generator.assert_called_once()
        read_parameters.assert_not_called()

    def test_missing_toppar_fails_before_fitting(self):
        (self.vanilla / "toppar").rename(self.vanilla / "saved_toppar")
        with self.assertRaisesRegex(FileNotFoundError, "toppar"):
            self.prepare()
        self.assertFalse(self.abfe.exists())
        self.generator.assert_not_called()

    def test_invalid_temperature_does_not_create_abfe_directory(self):
        for temperature in (0, -1, float("inf"), float("nan")):
            with self.subTest(temperature=temperature):
                with self.assertRaisesRegex(ValueError, "Temperature"):
                    self.prepare(temperature=temperature)
                self.assertFalse(self.abfe.exists())
                self.generator.assert_not_called()

    def test_analysis_failure_preserves_complex_and_working_directory(self):
        working_directory = Path.cwd()
        self.generator.side_effect = RuntimeError("analysis failed")

        with self.assertRaisesRegex(RuntimeError, "analysis failed"):
            self.prepare()
        self.assertTrue((self.abfe / "complex").is_dir())
        self.assertEqual(Path.cwd(), working_directory)


class BoreschTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.topology = self.root / "prod.tpr"
        self.trajectory = self.root / "prod.xtc"
        self.topology.touch()
        self.trajectory.touch()
        self.output = self.root / "boresch"

        self.universe = Mock()
        self.universe.select_atoms.return_value = SimpleNamespace(n_atoms=3)
        self.load_universe = self.enterContext(
            patch("vaasu.boresch.mda.Universe", return_value=self.universe)
        )
        self.find_ligand = self.enterContext(
            patch("vaasu.boresch.search.find_ligand_atoms", return_value=[[1, 2, 3]])
        )
        self.find_host = self.enterContext(patch("vaasu.boresch.search.FindHostAtoms"))
        self.find_host.return_value.host_atoms = [[4, 5, 6]]
        self.find_restraint = self.enterContext(
            patch("vaasu.boresch.FindBoreschRestraint")
        )
        self.restraint = self.find_restraint.return_value.restraint
        self.restraint.bond = SimpleNamespace(atype="bond", fc=10.0)
        self.restraint.angles = [
            SimpleNamespace(atype="angle", fc=10.0) for _ in range(2)
        ]
        self.restraint.dihedrals = [
            SimpleNamespace(atype="dihedral", fc=10.0) for _ in range(3)
        ]

    def generate(self, **kwargs):
        return generate_boresch_restraints(
            self.topology,
            self.trajectory,
            self.output,
            ligand_selection="resname unk",
            temperature=kwargs.get("temperature", 298.15),
        )

    def test_fits_before_writing_and_closes_trajectory(self):
        events = Mock()
        events.attach_mock(self.restraint.fit_fc, "fit")
        events.attach_mock(self.restraint.write, "write")
        events.attach_mock(self.restraint.plot, "plot")
        working_directory = Path.cwd()

        self.assertEqual(self.generate(), self.output)
        self.assertEqual(Path.cwd(), working_directory)
        self.find_ligand.assert_called_once_with(
            self.universe,
            l_selection="(resname unk) and not name H*",
            p_align="protein and name CA",
        )
        self.find_restraint.assert_called_once_with(
            self.universe, [([1, 2, 3], [4, 5, 6])]
        )
        self.assertEqual(
            events.mock_calls,
            [
                call.fit(temperature=298.15),
                call.write(path=str(self.output)),
                call.plot(path=str(self.output)),
            ],
        )
        self.universe.trajectory.close.assert_called_once()

    def test_existing_output_is_preserved_without_loading_trajectory(self):
        self.output.mkdir()
        sentinel = self.output / "keep.txt"
        sentinel.write_text("existing data")
        with self.assertRaises(FileExistsError):
            self.generate()
        self.assertEqual(sentinel.read_text(), "existing data")
        self.load_universe.assert_not_called()

    def test_invalid_inputs_do_not_load_trajectory_or_create_outputs(self):
        self.topology.unlink()
        with self.assertRaises(FileNotFoundError):
            self.generate()
        self.topology.touch()
        with self.assertRaisesRegex(ValueError, "Temperature"):
            self.generate(temperature=float("nan"))
        self.assertFalse(self.output.exists())
        self.load_universe.assert_not_called()

    def test_empty_selections_fail_before_creating_outputs(self):
        for selection in ("ligand", "host"):
            with self.subTest(selection=selection):
                self.universe.select_atoms.side_effect = [
                    SimpleNamespace(n_atoms=0 if selection == "ligand" else 3),
                    SimpleNamespace(n_atoms=0),
                ]
                with self.assertRaisesRegex(ValueError, selection):
                    self.generate()
                self.assertFalse(self.output.exists())
                self.find_ligand.assert_not_called()
        self.assertEqual(self.universe.trajectory.close.call_count, 2)

    def test_empty_ligand_or_host_candidates_raise_clear_error(self):
        self.find_ligand.return_value = []
        with self.assertRaisesRegex(ValueError, "No suitable"):
            self.generate()
        self.output.rmdir()
        self.find_ligand.return_value = [[1, 2, 3]]
        self.find_host.return_value.host_atoms = []
        with self.assertRaisesRegex(ValueError, "No suitable"):
            self.generate()
        self.find_restraint.assert_not_called()

    def test_invalid_fitted_constants_are_not_written(self):
        for index, force_constant in enumerate((0, -1, float("nan"), float("inf"))):
            with self.subTest(force_constant=force_constant):
                self.output = self.root / f"boresch_{index}"
                self.restraint.bond.fc = force_constant
                with self.assertRaisesRegex(ValueError, "fitted force constant"):
                    self.generate()
                self.restraint.write.assert_not_called()
                self.restraint.plot.assert_not_called()

    def test_write_failure_preserves_partial_output_and_closes_trajectory(self):
        def fail_write(*, path):
            (Path(path) / "partial.top").write_text("partial output")
            raise OSError("write failed")

        working_directory = Path.cwd()
        self.restraint.write.side_effect = fail_write
        with self.assertRaisesRegex(OSError, "write failed"):
            self.generate()
        self.assertEqual((self.output / "partial.top").read_text(), "partial output")
        self.assertEqual(Path.cwd(), working_directory)
        self.universe.trajectory.close.assert_called_once()
        self.restraint.plot.assert_not_called()


if __name__ == "__main__":
    unittest.main()
