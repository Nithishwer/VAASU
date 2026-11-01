"""Prepare one co-annihilation system for the complete lambda schedule."""

import re
from itertools import combinations_with_replacement
from math import fsum, sqrt
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

_INCLUDE = re.compile(r'^\s*#\s*include\s+["<]([^">]+)[">]', re.MULTILINE)


def _records(text: str):
    """Read section names and data without changing the original text."""
    section = ""
    for index, line in enumerate(text.splitlines()):
        content = line.partition(";")[0].strip()
        if not content:
            continue
        if content.startswith("["):
            section = content.strip("[] ").lower()
            yield index, section, []
        else:
            yield index, section, content.split()


def _read_parameters(
    root: Path, ligand_itp: str | Path, atomtype_itp: str | Path
) -> dict:
    """Read the supplied ITPs and build alchemical parameters in memory."""
    ligand_path = (root / ligand_itp).resolve()
    atomtype_path = (root / atomtype_itp).resolve()
    ligand_lines = ligand_path.read_text(encoding="utf-8").splitlines(keepends=True)
    atomtype_text = atomtype_path.read_text(encoding="utf-8")
    atoms = []
    for index, section, fields in _records("".join(ligand_lines)):
        if section == "atoms" and fields:
            atoms.append(fields)
            content, separator, comment = (
                ligand_lines[index].rstrip("\r\n").partition(";")
            )
            ligand_lines[index] = (
                f"{content.rstrip()} d{fields[1]} 0.0 {fields[7]}"
                + (f" ;{comment}" if separator else "")
                + "\n"
            )
    total_charge = fsum(float(atom[6]) for atom in atoms)
    charge = round(total_charge)
    if abs(total_charge - charge) > 1e-4:
        raise ValueError("Ligand charge must be within 1e-4 e of an integer.")
    atomtypes = {
        fields[0]: fields
        for _, section, fields in _records(atomtype_text)
        if section == "atomtypes" and fields
    }

    # Real-real and dummy-dummy overrides retain the legacy intraligand LJ terms.
    generated_types = atomtype_text.rstrip() + "\n"
    for fields in atomtypes.values():
        generated_types += (
            " ".join([f"d{fields[0]}", *fields[1:5], "0.0", "0.0"]) + "\n"
        )
    generated_types += "\n[ nonbond_params ]\n; i j funct sigma epsilon\n"
    for left, right in combinations_with_replacement(atomtypes.values(), 2):
        sigma = (float(left[5]) + float(right[5])) / 2
        epsilon = sqrt(float(left[6]) * float(right[6]))
        for prefix in ("", "d"):
            generated_types += (
                f"{prefix}{left[0]} {prefix}{right[0]} 1 {sigma:.8f} {epsilon:.8f}\n"
            )
    topology = (root / "complex.top").read_text(encoding="utf-8")
    molecules = {}
    for index, section, fields in _records(topology):
        if section == "molecules" and fields:
            molecules.setdefault(fields[0], []).append((index, int(fields[1])))
    return {
        "ligand_path": ligand_path,
        "atomtype_path": atomtype_path,
        "ligand_text": "".join(ligand_lines),
        "atomtype_text": generated_types,
        "charge": charge,
        "topology": topology,
        "molecules": molecules,
    }


def _prepare_counterions(frame_path: Path, ligand_selection: str, parameters: dict):
    """Convert existing counterions and return coordinates plus topology metadata."""
    text = frame_path.read_text(encoding="utf-8")
    charge = parameters["charge"]
    if charge == 0:
        return text, None

    lines = text.splitlines(keepends=True)
    universe = mda.Universe(str(frame_path), to_guess=())
    try:
        ligand = universe.select_atoms(ligand_selection)
        ordinary_name, atomtype, mass, ion_charge = (
            ("CL", "Cl", "35.45", -1) if charge > 0 else ("NA", "Na", "22.99", 1)
        )
        candidates = universe.select_atoms(f"resname {ordinary_name}")
        count = abs(charge)
        if candidates.n_atoms < count:
            raise ValueError(
                f"Need {count} {ordinary_name} counterions, found {candidates.n_atoms}."
            )
        distances = distance_array(ligand.center_of_geometry(), candidates.positions)[0]
        ranked = sorted(
            zip(distances, candidates.indices, strict=False),
            key=lambda item: (-item[0], item[1]),
        )
        selected = [int(index) for _, index in ranked[:count]]
        positions = [
            (0, 0, 0),
            (universe.dimensions[0] / 2, 0, 0),
            (0, universe.dimensions[1] / 2, 0),
        ]
        if count > len(positions):
            raise ValueError(
                "Alchemical-ion placement supports ligand charges from -3 to +3."
            )
        alchemical_name = f"B_{ordinary_name}"
        records = lines[2:-1]
        for index, position in zip(selected, positions, strict=False):
            record = records[index]
            decimals = [match.start() for match in re.finditer(r"\.", record[20:])]
            width = decimals[1] - decimals[0]
            # MDAnalysis uses angstroms; GRO positions are in nanometres.
            coordinates = "".join(
                f"{value / 10:{width}.{width - 5}f}" for value in position
            )
            records[index] = (
                record[:5]
                + f"{alchemical_name:<5}"
                + record[10:20]
                + coordinates
                + record[20 + 3 * width :]
            )
        selected_set = set(selected)
        order = [
            index for index in range(len(records)) if index not in selected_set
        ] + sorted(selected)
        coordinates = "".join(
            lines[:2] + [records[index] for index in order] + lines[-1:]
        )
        ion = {
            "ordinary_name": ordinary_name,
            "name": alchemical_name,
            "atomtype": atomtype,
            "mass": mass,
            "charge": ion_charge,
            "count": count,
            "remaining": candidates.n_atoms - count,
        }
        return coordinates, ion
    finally:
        universe.trajectory.close()


def _alchemical_path(path: Path) -> Path:
    return path.with_stem(f"{path.stem}_alchemical")


def _build_topology(root: Path, parameters: dict, ion: dict | None) -> str:
    topology = parameters["topology"]
    # Match include targets, never arbitrary occurrences in comments or other paths.
    for match in reversed(list(_INCLUDE.finditer(topology))):
        source = (root / match[1]).resolve()
        if source in (parameters["ligand_path"], parameters["atomtype_path"]):
            target = _alchemical_path(source).relative_to(root).as_posix()
            topology = topology[: match.start(1)] + target + topology[match.end(1) :]
    lines = topology.splitlines(keepends=True)
    records = list(_records(topology))
    restraint_lines = [
        index
        for index, section, fields in records
        if section == "intermolecular_interactions" and not fields
    ]
    system_line = next(
        index
        for index, section, fields in records
        if section == "system" and not fields
    )
    molecule_line = next(
        index
        for index, section, fields in records
        if section == "molecules" and not fields
    )
    if (
        len(restraint_lines) != 1
        or not system_line < molecule_line < restraint_lines[0]
    ):
        raise ValueError(
            "Expected [ system ], [ molecules ], then one [ intermolecular_interactions ] section."
        )
    if any(
        not fields and molecule_line < index < restraint_lines[0]
        for index, _, fields in records
    ):
        raise ValueError(
            "The molecule list must be followed by intermolecular restraints."
        )
    if ion is None:
        return topology
    ion_line = parameters["molecules"][ion["ordinary_name"]][0][0]
    _, separator, comment = lines[ion_line].rstrip("\r\n").partition(";")
    lines[ion_line] = (
        f"{ion['ordinary_name']} {ion['remaining']}"
        + (f" ;{comment}" if separator else "")
        + "\n"
    )
    # Ligand and counterion lose their charges together through coul-lambdas.
    definition = (
        f"\n[ moleculetype ]\n{ion['name']} 1\n[ atoms ]\n"
        f"1 {ion['atomtype']} 1 {ion['ordinary_name']} {ion['ordinary_name']} 1 "
        f"{ion['charge']} {ion['mass']} {ion['atomtype']} 0 {ion['mass']}\n"
        "[ position_restraints ]\n1 1 1000 1000 1000\n\n"
    )
    restraint_start = restraint_lines[0]
    # Keep the comments introducing the fitted restraints with their original block.
    while (
        restraint_start > molecule_line
        and not lines[restraint_start - 1].partition(";")[0].strip()
    ):
        restraint_start -= 1
    lines[restraint_start] = f"{ion['name']} {ion['count']}\n" + lines[restraint_start]
    lines[system_line] = definition + lines[system_line]
    return "".join(lines)


def prepare_alchemical_system(
    complex_path: str | Path,
    *,
    ligand_selection: str,
    ligand_itp: str | Path,
    atomtype_itp: str | Path,
) -> Path:
    """Write shared A/B endpoints; the MDP lambda schedule controls each leg."""
    complex_path = Path(complex_path).expanduser().resolve()
    parameters = _read_parameters(complex_path, ligand_itp, atomtype_itp)
    coordinates, ion = _prepare_counterions(
        complex_path / "complex.gro", ligand_selection, parameters
    )
    outputs = {
        complex_path / "complex_alchemical.gro": coordinates,
        complex_path / "complex_alchemical.top": _build_topology(
            complex_path, parameters, ion
        ),
        _alchemical_path(parameters["ligand_path"]): parameters["ligand_text"],
        _alchemical_path(parameters["atomtype_path"]): parameters["atomtype_text"],
    }
    # Check every output before writing; exclusive creation also protects against races.
    for path in outputs:
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Alchemical output already exists: {path}")
    for path, content in outputs.items():
        with path.open("x", encoding="utf-8") as output:
            output.write(content)
    return complex_path
