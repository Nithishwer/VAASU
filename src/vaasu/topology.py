"""Copy a vanilla topology and regroup trailing ions without moving anchors."""

import re
import shutil
from math import isfinite
from pathlib import Path

_INCLUDE = re.compile(r'^\s*#\s*include\s+["<]([^">]+)[">]', re.MULTILINE)


def validate_topology_inputs(vanilla_path: Path) -> list[Path]:
    """Check topology resources and return local files and directories to copy."""
    # Require the vanilla topology and its supporting parameter directory.
    topology_path = vanilla_path / "topol.top"
    if not topology_path.is_file():
        raise FileNotFoundError(f"Missing topology input: {topology_path}")
    if not (vanilla_path / "toppar").is_dir():
        raise FileNotFoundError(
            f"Missing topology directory: {vanilla_path / 'toppar'}"
        )

    # Validate molecule ordering before any expensive fitting.
    topology = topology_path.read_bytes().decode("utf-8")
    _regroup_topology(topology)
    resources = [Path("toppar")]
    # Inspect direct includes only; GROMACS handles full include preprocessing.
    for match in _INCLUDE.finditer(topology):
        include = Path(match[1])
        if not (vanilla_path / include).is_file():
            raise FileNotFoundError(
                f"Missing topology include: {vanilla_path / include}"
            )
        # External resources stay at their source paths instead of being copied.
        if include.is_absolute() or ".." in include.parts:
            continue

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
    molecules_end = None
    ions_started = False

    # Parse a comment-free view while retaining original lines for the output.
    for line in topology.splitlines(keepends=True):
        content, separator, comment = line.partition(";")
        content = content.strip()
        header = re.fullmatch(r"\[\s*([^\]]+?)\s*\]", content)
        if header:
            next_section = header[1].strip().lower()
            if next_section == "molecules":
                if molecules_found:
                    raise ValueError(
                        "Topology must contain a single [ molecules ] section."
                    )
                molecules_found = True
            elif section == "molecules":
                # Existing intermolecular restraints follow the molecule list.
                if next_section != "intermolecular_interactions":
                    raise ValueError(
                        "The [ molecules ] section must be the final system section."
                    )
                molecules_end = len(output)
            section = next_section
        elif section == "molecules" and content:
            fields = content.split()
            name = fields[0]
            if name in ion_counts:
                # Replace repeated ion entries with totals, keeping their comments.
                ions_started = True
                ion_counts[name] += int(fields[1])
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
    # Keep trailing comments and existing restraint text together after the ions.
    if molecules_end is None:
        molecules_end = len(output)
    while molecules_end and not output[molecules_end - 1].partition(";")[0].strip():
        molecules_end -= 1
    grouped = "".join(output[:molecules_end])
    if not grouped.endswith("\n"):
        grouped += "\n"
    for name, count in ion_counts.items():
        if count:
            grouped += f"{name:<12} {count}\n"
    return grouped + "".join(output[molecules_end:]), ion_counts


def _regroup_coordinates(coordinates: str) -> tuple[str, dict[str, int]]:
    """Reorder only trailing monatomic ions, preserving the original GRO records."""
    lines = coordinates.splitlines(keepends=True)
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
    # Keep parent-relative includes pointing to their original source after copying.
    for match in reversed(list(_INCLUDE.finditer(topology))):
        include = Path(match[1])
        if ".." in include.parts:
            source = (vanilla_path / include).resolve()
            topology = (
                topology[: match.start(1)] + str(source) + topology[match.end(1) :]
            )
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
