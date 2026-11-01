"""Check global lambda states, rendered simulation settings, and window outputs."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from vaasu.windows import create_windows


def read_mdp(path):
    settings = {}
    for line in path.read_text().splitlines():
        content = line.partition(";")[0]
        if "=" in content:
            key, value = content.split("=", 1)
            settings[key.strip().replace("_", "-").lower()] = value.split()
    return settings


class WindowTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.small_schedule = {
            "vdw": [0.0, 1.0],
            "coul": [0.0, 1.0],
            "rest": [0.0, 1.0],
        }

    def test_default_schedule_layout_and_stage_settings(self):
        working_directory = Path.cwd()
        windows = create_windows(self.root, temperature=298.15)
        expected_names = [
            f"{phase}.{number:02d}"
            for phase, count in (("rest", 12), ("coul", 11), ("vdw", 21))
            for number in range(count)
        ]
        self.assertEqual(windows, [self.root / name for name in expected_names])
        self.assertEqual(
            (self.root / "simulations_list.txt").read_text(),
            "".join(f"{name}\n" for name in expected_names),
        )
        stages = ("enmin", "nvt", "npt_b", "npt_pr", "prod")
        self.assertEqual(len(list(self.root.glob("*/*/*.mdp"))), 220)
        for state_index, window in enumerate(windows):
            self.assertEqual({path.name for path in window.iterdir()}, set(stages))
            phase, number = window.name.split(".")
            for stage in stages:
                settings = read_mdp(window / stage / f"{phase}.{stage}.{number}.mdp")
                self.assertEqual(settings["init-lambda-state"], [str(state_index)])
                for key in ("bonded-lambdas", "coul-lambdas", "vdw-lambdas"):
                    self.assertEqual(len(settings[key]), 44)
                if stage == "enmin":
                    self.assertEqual(settings["integrator"], ["steep"])
                    self.assertNotIn("ref-t", settings)
                    self.assertNotIn("gen-temp", settings)
                else:
                    self.assertEqual(settings["tc-grps"], ["Solute", "Solvent"])
                    self.assertEqual(settings["comm-grps"], ["Solute", "Solvent"])
                    self.assertEqual(settings["ref-t"], ["298.15", "298.15"])
        production = read_mdp(self.root / "rest.00/prod/rest.prod.00.mdp")
        arrays = {
            key: [float(value) for value in production[key]]
            for key in ("bonded-lambdas", "coul-lambdas", "vdw-lambdas")
        }
        self.assertEqual(
            arrays["bonded-lambdas"][:12],
            [0.0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.35, 0.5, 0.75, 1.0],
        )
        self.assertEqual(arrays["bonded-lambdas"][12:], [1.0] * 32)
        self.assertEqual(arrays["coul-lambdas"][:12], [0.0] * 12)
        self.assertEqual(
            arrays["coul-lambdas"][12:23], [index / 10 for index in range(11)]
        )
        self.assertEqual(arrays["coul-lambdas"][23:], [1.0] * 21)
        self.assertEqual(arrays["vdw-lambdas"][:23], [0.0] * 23)
        self.assertEqual(
            arrays["vdw-lambdas"][23:], [index / 20 for index in range(21)]
        )
        # Identical lambda vectors share the same Hamiltonian at leg boundaries.
        vectors = list(zip(*arrays.values(), strict=True))
        self.assertEqual(vectors[11], vectors[12])
        self.assertEqual(vectors[22], vectors[23])
        # NPT equilibration matches production geometry and reference pressure.
        self.assertEqual(production["nsteps"], ["10000000"])
        for stage, pressure in (
            ("npt_b", "Berendsen"),
            ("npt_pr", "Parrinello-Rahman"),
            ("prod", "Parrinello-Rahman"),
        ):
            settings = read_mdp(self.root / f"rest.00/{stage}/rest.{stage}.00.mdp")
            self.assertEqual(settings["pcoupl"], [pressure])
            self.assertEqual(settings["pcoupltype"], ["semiisotropic"])
            self.assertEqual(settings["ref-p"], ["1.0", "1.0"])
            self.assertEqual(settings["compressibility"], ["4.5e-05", "4.5e-05"])
        self.assertEqual(Path.cwd(), working_directory)

    def test_custom_schedule_sets_global_indices_arrays_and_temperature(self):
        windows = create_windows(
            self.root, temperature=310.0, lambda_schedule=self.small_schedule
        )
        self.assertEqual(
            [path.name for path in windows],
            ["rest.00", "rest.01", "coul.00", "coul.01", "vdw.00", "vdw.01"],
        )
        self.assertEqual(len(list(self.root.glob("*/*/*.mdp"))), 30)
        for state_index, window in enumerate(windows):
            phase, number = window.name.split(".")
            settings = read_mdp(window / "prod" / f"{phase}.prod.{number}.mdp")
            self.assertEqual(settings["init-lambda-state"], [str(state_index)])
            self.assertEqual(settings["ref-t"], ["310.0", "310.0"])
            self.assertEqual(settings["gen-temp"], ["310.0"])
            self.assertEqual(
                [float(value) for value in settings["bonded-lambdas"]],
                [0, 1, 1, 1, 1, 1],
            )
            self.assertEqual(
                [float(value) for value in settings["coul-lambdas"]], [0, 0, 0, 1, 1, 1]
            )
            self.assertEqual(
                [float(value) for value in settings["vdw-lambdas"]], [0, 0, 0, 0, 0, 1]
            )
        self.assertEqual(
            (self.root / "simulations_list.txt").read_text().splitlines(),
            [window.name for window in windows],
        )

    def test_preserves_supplied_order_duplicates_and_schedule(self):
        schedule = {
            "rest": [0.0, 0.7, 0.7, 0.4, 1.0],
            "coul": [0.0, 1.0],
            "vdw": [0.0, 1.0],
        }
        original = {phase: values.copy() for phase, values in schedule.items()}
        windows = create_windows(
            self.root, temperature=298.15, lambda_schedule=schedule
        )
        self.assertEqual(len(windows), 9)
        values = read_mdp(self.root / "rest.00/enmin/rest.enmin.00.mdp")[
            "bonded-lambdas"
        ]
        self.assertEqual(
            [float(value) for value in values[:5]], [0.0, 0.7, 0.7, 0.4, 1.0]
        )
        self.assertEqual(schedule, original)

    def test_checks_late_window_and_list_conflicts_before_creating_anything(self):
        conflict = self.root / "vdw.20"
        conflict.mkdir()
        sentinel = conflict / "keep.txt"
        sentinel.write_text("previous simulation")
        with self.assertRaises(FileExistsError):
            create_windows(self.root, temperature=298.15)
        self.assertEqual(list(self.root.iterdir()), [conflict])
        self.assertEqual(sentinel.read_text(), "previous simulation")
        sentinel.unlink()
        conflict.rmdir()
        conflict = self.root / "simulations_list.txt"
        conflict.write_text("previous list")
        with self.assertRaises(FileExistsError):
            create_windows(self.root, temperature=298.15)
        self.assertEqual(list(self.root.iterdir()), [conflict])
        self.assertEqual(conflict.read_text(), "previous list")
        conflict.unlink()
        conflict.symlink_to("missing.txt")
        with self.assertRaises(FileExistsError):
            create_windows(self.root, temperature=298.15)
        self.assertTrue(conflict.is_symlink())
        self.assertFalse((self.root / "rest.00").exists())

    def test_io_failure_preserves_completed_stages(self):
        original_open = Path.open

        def fail_nvt(path, *args, **kwargs):
            if path.name == "rest.nvt.00.mdp":
                raise OSError("write failed")
            return original_open(path, *args, **kwargs)

        with (
            patch.object(Path, "open", fail_nvt),
            self.assertRaisesRegex(OSError, "write failed"),
        ):
            create_windows(
                self.root, temperature=298.15, lambda_schedule=self.small_schedule
            )
        self.assertTrue((self.root / "rest.00/enmin/rest.enmin.00.mdp").is_file())
        self.assertFalse((self.root / "simulations_list.txt").exists())


if __name__ == "__main__":
    unittest.main()
