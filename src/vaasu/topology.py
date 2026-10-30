"""Copy a vanilla topology and regroup trailing ions without moving anchors."""

import re
import shutil
from math import isfinite
from pathlib import Path


def validate_topology_inputs(vanilla_path: Path) -> list[Path]:
    """Check the supported local include layout and return resources to copy."""
    # Require the vanilla topology and its supporting parameter directory.
    topology_path = vanilla_path / "topol.top"
    if not topology_path.is_file():
        raise FileNotFoundError(f"Missing topology input: {topology_path}")
    if not (vanilla_path / "toppar").is_dir():
        raise FileNotFoundError(
            f"Missing topology directory: {vanilla_path / 'toppar'}"
        )

    # Validate molecule ordering and existing restraints before any expensive fitting.
    topology = topology_path.read_bytes().decode("utf-8")
    _regroup_topology(topology)
    resources = [Path("toppar")]
    # Inspect direct includes only; GROMACS handles full include preprocessing.
    for line in topology.splitlines():
        if not re.match(r"\s*#\s*include\b", line):
            continue
        match = re.match(r'\s*#\s*include\s+"([^"]+)"', line)
        if match is None:
            raise ValueError(
                f"Only quoted local topology includes are supported: {line}"
            )
        # Preserve relative include paths without allowing source or output conflicts.
        include = Path(match[1])
        if include.is_absolute() or ".." in include.parts:
            raise ValueError(
                f"Topology include must stay beneath the vanilla run: {include}"
            )
        if include in {
            Path("complex.gro"),
            Path("complex.top"),
            Path("complex_ungrouped.gro"),
            Path("complex_ungrouped.top"),
            Path("BoreschRestraint.top"),
        }:
            raise ValueError(
                f"Topology include conflicts with a setup output: {include}"
            )
        if not (vanilla_path / include).is_file():
            raise FileNotFoundError(
                f"Missing topology include: {vanilla_path / include}"
            )

        # Copy whole force fields so their nested includes remain available.
        resource = include
        for index, part in enumerate(include.parts):
            if part.endswith(".ff"):
                resource = Path(*include.parts[: index + 1])
                break
        # Avoid copying files already covered by toppar or a force-field directory.
        if not any(resource.is_relative_to(existing) for existing in resources):
            resources.append(resource)
    return resources


def _regroup_topology(topology: str) -> tuple[str, dict[str, int]]:
    """Combine trailing NA/CL molecule entries, retaining other lines and comments."""
    output = []
    ion_counts = {"NA": 0, "CL": 0}
    section = ""
    molecules_found = False
    ions_started = False

    # Parse a comment-free view while retaining original lines for the output.
    for line in topology.splitlines(keepends=True):
        content, separator, comment = line.partition(";")
        content = content.strip()
        header = re.fullmatch(r"\[\s*([^\]]+?)\s*\]", content)
        if header:
            # Restraints will be appended after one final molecule list.
            section = header[1].strip().lower()
            if section == "intermolecular_interactions":
                raise ValueError(
                    "Input topology already contains intermolecular restraints."
                )
            if molecules_found:
                raise ValueError(
                    "The [ molecules ] section must be the final topology section."
                )
            molecules_found = section == "molecules"
        elif section == "molecules" and content:
            fields = content.split()
            try:
                if len(fields) != 2:
                    raise ValueError
                count = int(fields[1])
                if count < 0:
                    raise ValueError
            except ValueError:
                raise ValueError(f"Malformed molecule entry: {line.strip()}") from None
            name = fields[0]
            if name in ion_counts:
                # Replace repeated ion entries with totals, keeping their comments.
                ions_started = True
                ion_counts[name] += count
                if separator:
                    output.append(";" + comment)
                continue
            if ions_started:
                raise ValueError(
                    "NA/CL molecule entries must follow all non-ion molecules."
                )
        output.append(line)

    if not molecules_found:
        raise ValueError("Topology is missing its [ molecules ] section.")
    # Match the coordinate ordering: non-ion molecules, then NA, then CL.
    grouped = "".join(output)
    if not grouped.endswith("\n"):
        grouped += "\n"
    for name, count in ion_counts.items():
        if count:
            grouped += f"{name:<12} {count}\n"
    return grouped, ion_counts


def _regroup_coordinates(coordinates: str) -> tuple[str, dict[str, int]]:
    """Reorder only trailing monatomic ions, preserving the original GRO records."""
    lines = coordinates.splitlines(keepends=True)
    # Check the declared atom count and the three- or nine-value GRO box.
    try:
        atom_count = int(lines[1].strip())
        if atom_count < 0 or len(lines) != atom_count + 3:
            raise ValueError
        box = [float(value) for value in lines[-1].split()]
        if len(box) not in (3, 9) or not all(isfinite(value) for value in box):
            raise ValueError
    except (IndexError, ValueError):
        raise ValueError("Malformed GRO header, atom count, or box.") from None

    other_atoms = []
    ions = {"NA": [], "CL": []}
    ions_started = False
    for line in lines[2:-1]:
        # Validate fixed-width atom fields without rewriting coordinates or numbering.
        try:
            if len(line.rstrip("\r\n")) < 44:
                raise ValueError
            int(line[:5])
            int(line[15:20])
            position = [float(line[start : start + 8]) for start in (20, 28, 36)]
            if (
                not line[5:10].strip()
                or not line[10:15].strip()
                or not all(isfinite(value) for value in position)
            ):
                raise ValueError
        except ValueError:
            raise ValueError(f"Malformed GRO atom record: {line.strip()}") from None
        name = line[5:10].strip()
        if name in ions:
            ions_started = True
            ions[name].append(line)
        else:
            if ions_started:
                # Moving this atom ahead of an ion would change its restraint index.
                raise ValueError(
                    "NA/CL atoms must follow all non-ion atoms to preserve anchor indices."
                )
            other_atoms.append(line)

    # Keep the header, box and within-species atom order exactly as supplied.
    grouped = lines[:2] + other_atoms + ions["NA"] + ions["CL"] + lines[-1:]
    return "".join(grouped), {name: len(records) for name, records in ions.items()}


def prepare_topology(
    vanilla_path: str | Path,
    boresch_path: str | Path,
    complex_path: str | Path,
) -> Path:
    """Populate an existing empty complex directory using fresh restraint results."""
    vanilla_path = Path(vanilla_path).expanduser().resolve()
    boresch_path = Path(boresch_path).expanduser().resolve()
    complex_path = Path(complex_path).expanduser()
    # Only populate the empty directory created by the orchestrator.
    if complex_path.is_symlink():
        raise FileExistsError(f"Output already exists: {complex_path}")
    complex_path = complex_path.resolve()
    if not complex_path.is_dir():
        raise NotADirectoryError(f"Expected an empty complex directory: {complex_path}")
    if any(complex_path.iterdir()):
        raise FileExistsError(f"Complex directory is not empty: {complex_path}")

    # Check both source resources and fitted outputs before preparing any files.
    resources = validate_topology_inputs(vanilla_path)
    for name in ("ClosestRestraintFrame.gro", "BoreschRestraint.top"):
        if not (boresch_path / name).is_file():
            raise FileNotFoundError(f"Missing restraint output: {boresch_path / name}")

    topology = (vanilla_path / "topol.top").read_bytes().decode("utf-8")
    frame_path = boresch_path / "ClosestRestraintFrame.gro"
    coordinates = frame_path.read_bytes().decode("utf-8")
    # Coordinates and topology must describe the same number of each ion species.
    grouped_topology, topology_ions = _regroup_topology(topology)
    grouped_coordinates, coordinate_ions = _regroup_coordinates(coordinates)
    if coordinate_ions != topology_ions:
        raise ValueError(
            f"Ion counts disagree between coordinates {coordinate_ions} and topology {topology_ions}."
        )
    restraints = (boresch_path / "BoreschRestraint.top").read_bytes().decode("utf-8")

    # All checks above finish before copying anything into the complex directory.
    for resource in resources:
        source = vanilla_path / resource
        destination = complex_path / resource
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.is_dir():
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)
    # Keep the original ordering for inspection alongside the standalone restraints.
    shutil.copy2(frame_path, complex_path / "complex_ungrouped.gro")
    shutil.copy2(
        boresch_path / "BoreschRestraint.top", complex_path / "BoreschRestraint.top"
    )
    (complex_path / "complex_ungrouped.top").write_text(
        topology + "\n" + restraints, encoding="utf-8"
    )
    # Write matching grouped files and attach the fitted restraint text once.
    (complex_path / "complex.gro").write_text(grouped_coordinates, encoding="utf-8")
    (complex_path / "complex.top").write_text(
        grouped_topology + "\n" + restraints, encoding="utf-8"
    )
    return complex_path
