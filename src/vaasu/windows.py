"""Create lambda windows and render the packaged simulation-stage templates."""

from importlib.resources import files
from pathlib import Path

_DEFAULT_LAMBDA_SCHEDULE = {
    "rest": [0.0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.35, 0.5, 0.75, 1.0],
    "coul": [index / 10 for index in range(11)],
    "vdw": [index / 20 for index in range(21)],
}
_STAGES = ("enmin", "nvt", "npt_b", "npt_pr", "prod")


def create_windows(
    complex_path: str | Path,
    *,
    temperature: float,
    lambda_schedule: dict[str, list[float]] | None = None,
) -> list[Path]:
    """Write fresh window directories, MDPs, and their ordered simulation list."""
    complex_path = Path(complex_path).expanduser().resolve()
    schedule = _DEFAULT_LAMBDA_SCHEDULE if lambda_schedule is None else lambda_schedule
    # Keep repeated phase endpoints: every supplied value represents a window.
    states = (
        [
            ("rest", number, (value, 0.0, 0.0))
            for number, value in enumerate(schedule["rest"])
        ]
        + [
            ("coul", number, (1.0, value, 0.0))
            for number, value in enumerate(schedule["coul"])
        ]
        + [
            ("vdw", number, (1.0, 1.0, value))
            for number, value in enumerate(schedule["vdw"])
        ]
    )
    windows = [complex_path / f"{phase}.{number:02d}" for phase, number, _ in states]
    simulation_list = complex_path / "simulations_list.txt"
    for path in [*windows, simulation_list]:
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Window output already exists: {path}")

    # All stages and windows share the same global lambda arrays.
    lambda_arrays = {
        f"{name}_lambdas": " ".join(str(lambdas[column]) for _, _, lambdas in states)
        for column, name in enumerate(("bonded", "coul", "vdw"))
    }
    template_dir = files("vaasu").joinpath("templates", "mdp")
    templates = {
        stage: template_dir.joinpath(f"{stage}.mdp").read_text(encoding="utf-8")
        for stage in _STAGES
    }
    for state_index, (phase, number, _) in enumerate(states):
        window_path = windows[state_index]
        window_path.mkdir()
        for stage, template in templates.items():
            stage_path = window_path / stage
            stage_path.mkdir()
            mdp_path = stage_path / f"{phase}.{stage}.{number:02d}.mdp"
            with mdp_path.open("x", encoding="utf-8") as mdp:
                mdp.write(
                    template.format(
                        state_index=state_index,
                        temperature=temperature,
                        **lambda_arrays,
                    )
                )
    with simulation_list.open("x", encoding="utf-8") as output:
        output.writelines(f"{window.name}\n" for window in windows)

    print("A/B are global endpoints; each leg changes only its lambda component.")
    print(
        "Expected total charge: 0 e throughout, provided the input system is neutral."
    )
    print(f"{'Window':<12} {'State':>5} {'Bonded':>8} {'Coul':>8} {'VDW':>8}")
    for state_index, (phase, number, lambdas) in enumerate(states):
        name = f"{phase}.{number:02d}"
        print(
            f"{name:<12} {state_index:>5}"
            + "".join(f"{value:>8g}" for value in lambdas)
        )
    return windows
