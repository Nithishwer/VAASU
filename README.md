# VAASU

VAASU's An ABFEP Setup Utility

## Setup and activation

```bash
cd /biggin/b211/reub0138/Projects/VAASU
uv sync --locked
source .venv/bin/activate
```

This creates `.venv` and installs VAASU in editable mode with Ruff.

MDRestraintsGenerator and its Python dependencies install automatically from the
pinned [fork](https://github.com/Nithishwer/MDRestraintsGenerator); no separate clone is needed.
Once activated, use `python` directly from any directory.

```bash
python -c "import sys; print(sys.executable)"  # Check the interpreter
python /absolute/path/to/script.py
deactivate
```

## Run without activation

From the repository:

```bash
uv run python
uv run python /absolute/path/to/script.py
```

From another directory:

```bash
uv run --project /biggin/b211/reub0138/Projects/VAASU python /absolute/path/to/script.py
```

## Prepare complex directories and restraints

Each ligand has an MDAnalysis selection and a list of vanilla runs. Each run
specifies its fitting temperature in kelvin and one or more ABFE destinations:

```yaml
ligands:
  - name: DOR_10
    ligand_selection: "resname unk"
    runs:
      - temperature: 298.15
        vanilla_path: /data/DOR_10/vanilla_rep_1
        abfe_paths:
          - /data/DOR_10/abfe_van1_hrex_r1
```

Use the configuration to prepare each destination:

```python
from vaasu.complex import prepare_complex
from vaasu.config import load_config

config = load_config("abfe_campaign.yaml")
for ligand in config["ligands"]:
    for run in ligand["runs"]:
        for abfe_path in run["abfe_paths"]:
            complex_path = prepare_complex(
                run["vanilla_path"],
                abfe_path,
                ligand_selection=ligand["ligand_selection"],
                temperature=run["temperature"],
            )
            print(complex_path)
```

The vanilla directory must contain compatible `prod.tpr` and `prod.xtc` files,
`topol.top`, and a `toppar/` directory. Quoted topology includes must use local
paths beneath the vanilla directory; absolute paths and `..` are rejected.
Directly referenced include files are checked before restraint fitting. Referenced
local `.ff` directories and the whole `toppar/` directory are copied, retaining
their nested includes. Other directly referenced files are copied at their
relative paths. Full include preprocessing remains GROMACS's responsibility.
Restraint searches use protein CA atoms as host anchors and exclude ligand atoms
whose names begin with H. No GROMACS executable or separate generator script is
needed for this step.

Each ABFE destination receives a `boresch/` directory containing
`ClosestRestraintFrame.gro`, `BoreschRestraint.top`, and diagnostic plots. Its sibling
`complex/` contains:

- `complex.gro`: the fitted frame with trailing ions grouped as NA, then CL.
- `complex.top`: the vanilla topology with matching aggregated ion counts and the
  generated Boresch restraint text appended once.
- `complex_ungrouped.gro` and `complex_ungrouped.top`: the original frame and
  topology order, with restraints attached to the topology.
- `BoreschRestraint.top`, `toppar/`, and the referenced force-field resources.

All NA/CL atoms and molecule entries must follow the non-ion system. Regrouping
preserves non-ion atom indices, atom records, coordinates, and the box; ion counts
must agree between coordinates and topology. Inputs already containing an
`[ intermolecular_interactions ]` section are rejected. The topology must have a
single final `[ molecules ]` section. Alchemical systems, indices, simulation
windows, and job scripts are not generated yet.

Existing `complex` or `boresch` paths, including symlinks, raise `FileExistsError`.
Failures propagate and partial output is preserved. Restraint failures leave the
empty `complex/` directory and any partial `boresch/` output. Inspect and remove or
move partial directories explicitly before retrying. Source files are not changed.

The Yuan example can be run from any working directory after activating VAASU:

```bash
python /biggin/b211/reub0138/Projects/orexin/yuan_DOR_c1_v2_w_na/s2_vaasu_setup_abfe.py
```

## Development

Run from the repository; replace `<package>` with the dependency name.

```bash
uv add <package>          # Runtime dependency
uv add --dev <package>    # Development dependency
uv remove <package>
uv run ruff check .
uv run ruff format .
uv run python -m unittest discover -s tests
uv build                 # Wheel and source distribution in dist/
```

Commit `pyproject.toml` and `uv.lock` together after dependency changes.
