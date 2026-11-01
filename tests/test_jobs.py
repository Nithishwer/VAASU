"""Check ARCHER2 resources, window order, and generated job execution flow."""

import json
import os
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from vaasu.jobs import write_job_scripts
from vaasu.windows import create_windows


class JobTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory(prefix="vaasu jobs ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.windows = create_windows(
            self.root,
            temperature=298.15,
            lambda_schedule={"rest": [0, 1], "coul": [0, 1], "vdw": [0, 1]},
        )
        for filename in (
            "complex_alchemical.gro",
            "complex_alchemical.top",
            "index_alchemical.ndx",
        ):
            (self.root / filename).write_text("input fixture\n")
        self.events_path = self.root / "events.jsonl"
        binary_dir = self.root / "bin"
        binary_dir.mkdir()
        # Emulate job steps locally: no Slurm or GROMACS executable is invoked.
        srun = binary_dir / "srun"
        srun.write_text("""#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

args = sys.argv[1:]
with Path(os.environ["JOB_EVENTS"]).open("a") as output:
    output.write(json.dumps({"cwd": str(Path.cwd()), "args": args, "omp": os.environ["OMP_NUM_THREADS"], "libraries": os.environ["LD_LIBRARY_PATH"]}) + "\\n")
operation = args[args.index("gmx_mpi") + 1]
if operation == "grompp":
    if os.environ.get("FAIL_GROMPP"):
        sys.exit(7)
    for option in ("-f", "-c", "-r", "-n", "-p"):
        assert Path(args[args.index(option) + 1]).is_file(), args
    Path(args[args.index("-o") + 1]).write_text("generated TPR\\n")
else:
    assert operation == "mdrun", args
    directories = []
    for value in args[args.index("-multidir") + 1:]:
        if value.startswith("-"):
            break
        directories.append(Path(value))
    for directory in directories:
        assert (directory / "topol.tpr").is_file(), args
        if "-cpi" in args:
            assert (directory / args[args.index("-cpi") + 1]).is_file(), args
        (directory / "confout.gro").write_text("stage coordinates\\n")
        (directory / "state.cpt").write_text("checkpoint\\n")
    if "-plumed" in args:
        assert Path(args[args.index("-plumed") + 1]).is_file(), args
""")
        srun.chmod(0o755)
        self.environment = {
            **os.environ,
            "PATH": str(binary_dir) + os.pathsep + os.environ["PATH"],
            "SLURM_SUBMIT_DIR": str(self.root),
            "SLURM_CPUS_PER_TASK": "64",
            "JOB_EVENTS": str(self.events_path),
        }
        self.environment.pop("LD_LIBRARY_PATH", None)

    def write_jobs(self, windows=None):
        return write_job_scripts(
            self.root,
            windows=self.windows if windows is None else windows,
            job_name="DOR_10",
        )

    def run_job(self, name, **environment):
        return subprocess.run(
            ["bash", str(self.root / name)],
            cwd=self.root.parent,
            check=False,
            env={**self.environment, **environment},
            capture_output=True,
            text=True,
        )

    def events(self):
        return [json.loads(line) for line in self.events_path.read_text().splitlines()]

    def test_writes_five_variants_and_reference_resources(self):
        windows = [
            self.root / f"{phase}.{number:02d}"
            for phase, count in (("rest", 12), ("coul", 11), ("vdw", 21))
            for number in range(count)
        ]
        working_directory = Path.cwd()
        paths = self.write_jobs(windows)
        expected = {
            "job_complex_archer_test.sh": ("short", "00:20:00", 0),
            "job_complex_archer.sh": ("standard", "24:00:00", 0),
            "job_complex_archer_contd.sh": ("standard", "24:00:00", 1),
            "job_complex_archer_lowpriority.sh": ("lowpriority", "24:00:00", 0),
            "job_complex_archer_contd_lowpriority.sh": ("lowpriority", "24:00:00", 1),
        }
        self.assertEqual(paths, [self.root / name for name in expected])
        for path in paths:
            script = path.read_text()
            qos, walltime, continuation = expected[path.name]
            for directive in (
                "--job-name=DOR_10",
                "--nodes=22",
                "--ntasks=44",
                "--ntasks-per-node=2",
                "--cpus-per-task=64",
                "--account=e280-biggin",
                "--partition=standard",
                f"--qos={qos}",
                f"--time={walltime}",
            ):
                self.assertIn(f"#SBATCH {directive}\n", script)
            self.assertIn(f"continuation={continuation}\n", script)
            self.assertIn(
                "windows=(" + " ".join(window.name for window in windows) + ")", script
            )
            self.assertNotRegex(script, r"@[a-z_]+@")
            syntax = subprocess.run(
                ["bash", "-n", str(path)], check=False, capture_output=True, text=True
            )
            self.assertEqual(syntax.returncode, 0, syntax.stderr)
        self.assertEqual((self.root / "plumed.dat").read_bytes(), b"")
        self.assertEqual(Path.cwd(), working_directory)

    def test_sizes_custom_and_odd_window_counts(self):
        for count, nodes, tasks_per_node, cpus in (
            (6, 3, 2, 64),
            (9, 5, 2, 64),
            (1, 1, 1, 128),
        ):
            with self.subTest(count=count), TemporaryDirectory() as temporary:
                root = Path(temporary)
                paths = write_job_scripts(
                    root,
                    windows=[root / f"rest.{number:02d}" for number in range(count)],
                    job_name="custom",
                )
                script = paths[0].read_text()
                for directive in (
                    f"--ntasks={count}",
                    f"--nodes={nodes}",
                    f"--ntasks-per-node={tasks_per_node}",
                    f"--cpus-per-task={cpus}",
                ):
                    self.assertIn(f"#SBATCH {directive}\n", script)

    def test_initial_job_runs_stages_and_paths_in_supplied_order(self):
        windows = list(reversed(self.windows))
        self.write_jobs(windows)
        result = self.run_job("job_complex_archer_test.sh")
        self.assertEqual(result.returncode, 0, result.stderr)
        events = self.events()
        self.assertEqual(len(events), 35)  # Six grompp calls and one mdrun per stage.
        previous = None
        for stage_number, stage in enumerate(
            ("enmin", "nvt", "npt_b", "npt_pr", "prod")
        ):
            preprocessing = events[stage_number * 7 : stage_number * 7 + 6]
            for event, window in zip(preprocessing, windows, strict=True):
                args = event["args"]
                phase, number = window.name.split(".")
                self.assertEqual(event["cwd"], str(window / stage))
                self.assertEqual(event["omp"], "1")
                self.assertIn("--ntasks=1", args)
                self.assertIn("--cpus-per-task=1", args)
                self.assertEqual(
                    args[args.index("-f") + 1], f"{phase}.{stage}.{number}.mdp"
                )
                expected_coordinates = (
                    self.root / "complex_alchemical.gro"
                    if previous is None
                    else window / previous / "confout.gro"
                )
                self.assertEqual(args[args.index("-c") + 1], str(expected_coordinates))
                self.assertEqual(args[args.index("-r") + 1], str(expected_coordinates))
                self.assertEqual(
                    args[args.index("-n") + 1], str(self.root / "index_alchemical.ndx")
                )
                self.assertEqual(
                    args[args.index("-p") + 1],
                    str(self.root / "complex_alchemical.top"),
                )
            event = events[stage_number * 7 + 6]
            args = event["args"]
            self.assertEqual(event["omp"], "64")
            self.assertIn("--ntasks=6", args)
            start = args.index("-multidir") + 1
            self.assertEqual(
                args[start : start + 6],
                [f"{window.name}/{stage}" for window in windows],
            )
            self.assertNotIn("-cpi", args)
            if stage == "prod":
                self.assertEqual(
                    args[-7:],
                    [
                        "-replex",
                        "200",
                        "-hrex",
                        "-plumed",
                        str(self.root / "plumed.dat"),
                        "-dlb",
                        "no",
                    ],
                )
            else:
                self.assertNotIn("-hrex", args)
            self.assertIn("gromacs-2020.6/2020.6+plumed/lib", event["libraries"])
            previous = stage

    def test_preserves_three_digit_window_numbers(self):
        window = self.root / "rest.100"
        for stage in ("enmin", "nvt", "npt_b", "npt_pr", "prod"):
            directory = window / stage
            directory.mkdir(parents=True)
            (directory / f"rest.{stage}.100.mdp").write_text("MDP fixture\n")
        self.write_jobs([window, self.windows[2]])
        result = self.run_job("job_complex_archer.sh")
        self.assertEqual(result.returncode, 0, result.stderr)
        filenames = [
            event["args"][event["args"].index("-f") + 1]
            for event in self.events()
            if "grompp" in event["args"] and Path(event["cwd"]).parent == window
        ]
        self.assertEqual(
            filenames,
            [
                f"rest.{stage}.100.mdp"
                for stage in ("enmin", "nvt", "npt_b", "npt_pr", "prod")
            ],
        )

    def test_continuation_uses_only_production_checkpoints(self):
        self.write_jobs()
        result = self.run_job("job_complex_archer.sh")
        self.assertEqual(result.returncode, 0, result.stderr)
        before = len(self.events())
        result = self.run_job("job_complex_archer_contd.sh")
        self.assertEqual(result.returncode, 0, result.stderr)
        events = self.events()[before:]
        self.assertEqual(len(events), 1)
        args = events[0]["args"]
        self.assertIn("mdrun", args)
        self.assertNotIn("grompp", args)
        self.assertEqual(args[-2:], ["-cpi", "state.cpt"])
        start = args.index("-multidir") + 1
        self.assertEqual(
            args[start : start + 6], [f"{window.name}/prod" for window in self.windows]
        )

    def test_preprocessing_failure_stops_the_job(self):
        self.write_jobs()
        result = self.run_job("job_complex_archer.sh", FAIL_GROMPP="1")
        self.assertEqual(result.returncode, 7)
        self.assertEqual(len(self.events()), 1)
        self.assertNotIn("mdrun", self.events()[0]["args"])

    def test_refuses_conflicting_scripts_and_preserves_plumed(self):
        conflict = self.root / "job_complex_archer_contd_lowpriority.sh"
        conflict.write_text("existing job")
        with self.assertRaises(FileExistsError):
            self.write_jobs()
        self.assertEqual(conflict.read_text(), "existing job")
        self.assertFalse((self.root / "job_complex_archer_test.sh").exists())
        self.assertFalse((self.root / "plumed.dat").exists())
        conflict.unlink()
        conflict.symlink_to("missing.sh")
        with self.assertRaises(FileExistsError):
            self.write_jobs()
        self.assertTrue(conflict.is_symlink())
        conflict.unlink()
        (self.root / "plumed.dat").write_text("PRINT ARG=* FILE=colvar\n")
        self.write_jobs()
        self.assertEqual(
            (self.root / "plumed.dat").read_text(), "PRINT ARG=* FILE=colvar\n"
        )


if __name__ == "__main__":
    unittest.main()
