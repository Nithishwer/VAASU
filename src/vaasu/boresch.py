"""Generate Boresch restraints from a production trajectory."""

from math import isfinite
from pathlib import Path

import MDAnalysis as mda
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindBoreschRestraint


def generate_boresch_restraints(
    topology: str | Path,
    trajectory: str | Path,
    output_dir: str | Path,
    *,
    ligand_selection: str,
    temperature: float,
    host_selection: str = "protein and name CA",
) -> Path:
    """Write fitted restraints and diagnostic plots to a fresh directory."""
    topology = Path(topology).expanduser().resolve()
    trajectory = Path(trajectory).expanduser().resolve()
    output_dir = Path(output_dir).expanduser()

    # Check inputs and reject existing outputs before loading the trajectory.
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f"Restraint output already exists: {output_dir}")
    output_dir = output_dir.resolve()
    for input_path in (topology, trajectory):
        if not input_path.is_file():
            raise FileNotFoundError(f"Missing production input: {input_path}")
    if not isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be finite and positive, in kelvin.")

    universe = mda.Universe(str(topology), str(trajectory))
    try:
        # Use heavy ligand atoms and require nonempty ligand and host selections.
        heavy_ligand_selection = f"({ligand_selection}) and not name H*"
        if not universe.select_atoms(heavy_ligand_selection).n_atoms:
            raise ValueError(
                f"No heavy atoms match ligand selection: {ligand_selection}"
            )
        if not universe.select_atoms(host_selection).n_atoms:
            raise ValueError(f"No atoms match host selection: {host_selection}")

        output_dir.mkdir(parents=True, exist_ok=False)
        # Find ligand anchor triplets, then suitable host anchors for each candidate.
        ligand_candidates = search.find_ligand_atoms(
            universe,
            l_selection=heavy_ligand_selection,
            p_align=host_selection,
        )
        candidate_pairs = []
        for ligand_atoms in ligand_candidates:
            host_search = search.FindHostAtoms(
                universe, ligand_atoms[0], p_selection=host_selection
            )
            host_search.run()
            candidate_pairs.extend(
                (ligand_atoms, host_atoms) for host_atoms in host_search.host_atoms
            )
        if not candidate_pairs:
            raise ValueError("No suitable Boresch restraint anchors found.")

        # Evaluate candidates over the trajectory and fit at the requested temperature.
        analysis = FindBoreschRestraint(universe, candidate_pairs)
        analysis.run()
        restraint = analysis.restraint
        restraint.fit_fc(temperature=temperature)
        # Reject invalid fitted constants before writing restraints or diagnostic plots.
        for component in (restraint.bond, *restraint.angles, *restraint.dihedrals):
            if not isfinite(component.fc) or component.fc <= 0:
                raise ValueError(
                    f"Invalid fitted force constant for {component.atype}: {component.fc}"
                )
        restraint.write(path=str(output_dir))
        restraint.plot(path=str(output_dir))
    finally:
        # Release the trajectory even when fitting or writing fails.
        universe.trajectory.close()

    return output_dir
