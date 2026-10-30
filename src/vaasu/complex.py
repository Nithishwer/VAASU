"""Prepare a complex directory and its Boresch restraints."""

from math import isfinite
from pathlib import Path

from vaasu.boresch import generate_boresch_restraints
from vaasu.topology import prepare_topology, validate_topology_inputs


def validate_inputs(
    vanilla_path: Path,
    abfe_path: Path,
    *,
    temperature: float,
) -> None:
    """Check setup inputs and output conflicts before creating any directories."""
    # Reject existing outputs to protect earlier runs.
    for name in ("complex", "boresch"):
        output_path = abfe_path / name
        if output_path.exists():
            raise FileExistsError(f"Output already exists: {output_path}")

    # Check all required inputs before starting the restraint search.
    for name in ("prod.tpr", "prod.xtc"):
        input_path = vanilla_path / name
        if not input_path.is_file():
            raise FileNotFoundError(f"Missing production input: {input_path}")
    if not isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be finite and positive, in kelvin.")
    validate_topology_inputs(vanilla_path)


def prepare_complex(
    vanilla_path: str | Path,
    abfe_path: str | Path,
    *,
    ligand_selection: str,
    temperature: float,
) -> Path:
    """Generate restraints and prepare the complex coordinates and topology."""
    # Resolve paths once so every step works independently of the working directory.
    vanilla_path = Path(vanilla_path).expanduser().resolve()
    abfe_path = Path(abfe_path).expanduser().resolve()
    complex_path = abfe_path / "complex"
    boresch_path = abfe_path / "boresch"

    validate_inputs(vanilla_path, abfe_path, temperature=temperature)

    # Keep the complex directory and partial restraints available if fitting fails.
    complex_path.mkdir(parents=True, exist_ok=False)
    generate_boresch_restraints(
        vanilla_path / "prod.tpr",
        vanilla_path / "prod.xtc",
        boresch_path,
        ligand_selection=ligand_selection,
        temperature=temperature,
    )
    # Build the topology using the fitted frame, which matches the restraint indices.
    return prepare_topology(vanilla_path, boresch_path, complex_path)
