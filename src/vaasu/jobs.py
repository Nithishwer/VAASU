"""Write the reference ARCHER2 HREX jobs using the prepared window paths."""

from importlib.resources import files
from pathlib import Path
from shlex import join, quote

_JOBS = (
    ("job_complex_archer_test.sh", "short", "00:20:00", False),
    ("job_complex_archer.sh", "standard", "24:00:00", False),
    ("job_complex_archer_contd.sh", "standard", "24:00:00", True),
    ("job_complex_archer_lowpriority.sh", "lowpriority", "24:00:00", False),
    ("job_complex_archer_contd_lowpriority.sh", "lowpriority", "24:00:00", True),
)


def write_job_scripts(
    complex_path: str | Path,
    *,
    windows: list[Path],
    job_name: str,
) -> list[Path]:
    """Create five fresh job scripts and an empty PLUMED file if absent."""
    complex_path = Path(complex_path).expanduser().resolve()
    window_names = [window.relative_to(complex_path).as_posix() for window in windows]
    # Match the reference's two replicas per 128-core node, with one MPI rank each.
    tasks_per_node = min(2, len(windows))
    replacements = {
        "job_name": quote(job_name),
        "nodes": (len(windows) + tasks_per_node - 1) // tasks_per_node,
        "ntasks": len(windows),
        "tasks_per_node": tasks_per_node,
        "cpus_per_task": 128 // tasks_per_node,
        "windows": join(window_names),
    }
    template = (
        files("vaasu")
        .joinpath("templates", "jobs", "archer.sh")
        .read_text(encoding="utf-8")
    )
    outputs = {}
    for filename, qos, time_limit, continuation in _JOBS:
        content = template
        values = {
            **replacements,
            "qos": qos,
            "time_limit": time_limit,
            "continuation": int(continuation),
        }
        # Distinct markers leave Bash variables and array expressions untouched.
        for name, value in values.items():
            content = content.replace(f"@{name}@", str(value))
        outputs[complex_path / filename] = content
    for path in outputs:
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Job script already exists: {path}")
    plumed_path = complex_path / "plumed.dat"
    if not plumed_path.exists():
        with plumed_path.open("x", encoding="utf-8"):
            pass
    for path, content in outputs.items():
        with path.open("x", encoding="utf-8") as output:
            output.write(content)
    return list(outputs)
