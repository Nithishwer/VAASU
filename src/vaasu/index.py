"""Write shared coupling groups from the alchemical coordinates."""

from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.selections.gromacs import SelectionWriter


def write_index_file(
    complex_path: str | Path,
    *,
    solvent_selection: str = "resname SOL NA CL B_NA B_CL",
) -> Path:
    """Write Solute and Solvent groups for every lambda window."""
    complex_path = Path(complex_path).expanduser().resolve()
    path = complex_path / "index_alchemical.ndx"
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Index output already exists: {path}")

    universe = mda.Universe(str(complex_path / "complex_alchemical.gro"), to_guess=())
    try:
        solvent = universe.select_atoms(solvent_selection)
        solute = universe.select_atoms(f"not ({solvent_selection})")
        with SelectionWriter(str(path)) as writer:
            writer.write(solute, name="Solute")
            writer.write(solvent, name="Solvent")
    finally:
        universe.trajectory.close()
    return path
