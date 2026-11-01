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

## Prepare complex alchemical systems

Each ligand has an MDAnalysis selection and a list of vanilla runs. Each run
specifies its fitting temperature in kelvin and one or more ABFE destinations:

```yaml
ligands:
  - name: DOR_10
    ligand_selection: "resname unk"
    ligand_itp: toppar/10.itp
    atomtype_itp: toppar/10_atomtypes.itp
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
                ligand_itp=ligand["ligand_itp"],
                atomtype_itp=ligand["atomtype_itp"],
                solvent_selection=config.get(
                    "solvent_selection", "resname SOL NA CL B_NA B_CL"
                ),
                lambda_schedule=config.get("lambda_schedule"),
                job_name=ligand["name"],
            )
            print(complex_path)
```

The vanilla directory must contain compatible `prod.tpr` and `prod.xtc` files,
`topol.top`, and a `toppar/` directory. Includes can use quotes or angle brackets.
Directly referenced include files are checked before restraint fitting. Referenced
local `.ff` directories and the whole `toppar/` directory are copied, retaining
their nested includes. Other directly referenced files are copied at their
relative paths. Absolute includes remain external; includes containing `..` are
resolved to absolute source paths in the copied topology. Full include
preprocessing remains GROMACS's responsibility.
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
must agree between coordinates and topology. The topology must have a single
`[ molecules ]` section as its final system section; existing intermolecular
restraint sections are preserved after it. Generated restraints are appended.
The setup also writes one shared `complex_alchemical.gro/.top`, with
`*_alchemical.itp` ligand and atomtype files beside the originals. It then writes
`index_alchemical.ndx`, lambda-window directories and MDPs,
`simulations_list.txt`, and five ARCHER2 job scripts.

Both ITP paths are required keyword arguments to `prepare_complex()`. They are
relative to the vanilla directory and retain the same relative location in the
prepared complex. Existing callers and campaign files must add these two values.
The external Yuan example shown below must also be updated before running it.

### Alchemical protocol

VAASU implements co-annihilation with one topology for the entire lambda schedule.
The generator assumes Amber-style combination rule 2, one ligand residue with
unperturbed eight-column atom records, and seven-column atomtype definitions.
Supply matching ITPs referenced by the topology. Parameter files are read during
alchemical generation; there is no separate alchemical preflight validation.
Full GROMACS include preprocessing remains outside VAASU.

Ligand charge is the rounded sum of the parameterized atom charges, which must
be within `1e-4 e` of an integer. Each ligand atom retains its A-state parameters
and receives a B-state dummy type with zero
charge and unchanged mass. Dummy types have zero LJ parameters; real-real and
dummy-dummy nonbonded overrides retain the legacy intraligand LJ terms.

**A and B describe global endpoints, not the start and end of each leg.** The MDP
lambda components determine which parameters are active. All windows use the same
topology, coordinates, and index groups as their setup inputs; later simulation
stages use each window's preceding coordinates.

For charged ligands, VAASU converts the required number of opposite-charge NA/CL
ions to B_NA/B_CL. It selects the farthest ions from the ligand center of geometry
using non-periodic distances, with original atom order resolving ties. The legacy
placement positions are `(0, 0, 0)`, `(Lx/2, 0, 0)`, and `(0, Ly/2, 0)`; selected
ions move to the end of the coordinates. Protein and ligand atom ordering, other
coordinate records, box text, and fitted restraints are preserved. Alchemical ions
retain their LJ type and receive the legacy position restraints with force
constants of `1000 kJ mol^-1 nm^-2`. Placement supports ligand charges from -3 to
+3 and requires enough existing counterions.

| Parameter | Global A | Global B |
| --- | --- | --- |
| Ligand charges | Original | Zero |
| Alchemical-ion charges | Original | Zero |
| Ligand LJ types | Original | Dummy |
| Alchemical-ion LJ types | Original | Original |
| Boresch restraint force constants | Zero | Fitted |

During `rest`, Coulomb and VDW lambdas remain zero, so the ligand and counterions
retain their original charges and LJ interactions while the Boresch restraints
turn on. During `coul`, ligand and counterion charges turn off together. During
`vdw`, both stay uncharged while ligand LJ interactions turn off. Total charge
therefore stays equal to the input system's charge, within charge-rounding
tolerance; a neutral input remains neutral.

Neutral ligands need no alchemical ions and retain their original coordinates.
Counterion preparation uses the supplied coordinates and ligand selection directly.
The code does not infer full-system charge from residue names. Check full-system
A/B charges and preprocessing with GROMACS before simulation.

To add the alchemical system to an already prepared complex, call this step directly:

```python
from vaasu.alchemy import prepare_alchemical_system

prepare_alchemical_system(
    "/data/DOR_10/abfe_van1_hrex_r1/complex",
    ligand_selection="resname unk",
    ligand_itp="toppar/10.itp",
    atomtype_itp="toppar/10_atomtypes.itp",
)
```

Alchemical outputs are built in memory before writing and are created exclusively.
Existing output files or symlinks raise `FileExistsError`. Input and source files
are preserved; an I/O failure leaves partial output available for inspection.

This replaces the old phase-specific outputs and the standalone helpers are now
`prepare_alchemical_system()` and `write_index_file()`. Regenerate the system,
index, job scripts, and TPRs together in a fresh destination when updating an old
setup. Old production TPRs and continuation checkpoints retain the previous
topologies.

### Index files and lambda windows

After alchemical preparation, `prepare_complex()` writes `index_alchemical.ndx`
from `complex_alchemical.gro`. Every window uses these same two groups, with
GROMACS one-based atom numbers:

- `Solvent`: atoms matching `resname SOL NA CL B_NA B_CL` by default.
- `Solute`: every other atom, including protein, ligand, and membrane.

The optional `solvent_selection` argument accepts an MDAnalysis selection. Both
groups are referenced by name in the MDPs. No GROMACS executable is required.

The default lambda schedule reproduces the legacy 44 windows:

| Phase | Window directories | Global state indices |
| --- | --- | --- |
| rest | `rest.00` through `rest.11` | 0–11 |
| coul | `coul.00` through `coul.10` | 12–22 |
| vdw | `vdw.00` through `vdw.20` | 23–43 |

Each window has `enmin`, `nvt`, `npt_b`, `npt_pr`, and `prod` directories, with MDP
names such as `rest.00/nvt/rest.nvt.00.mdp`. The default produces 220 MDP files.
Five stage templates ship with VAASU. Both NPT equilibration stages use
semi-isotropic pressure coupling at 1.0 bar in xy and z, matching production's
coupling geometry and reference pressure. The legacy barostats and durations
remain unchanged: 10 ps NVT, 100 ps restrained NPT, and 500 ps unrestrained NPT.
These durations do not establish sufficient membrane or alchemical-window
equilibration. Template changes apply to newly generated windows; existing MDPs
must be updated and TPRs regenerated before rerunning the affected stages.
Every MDP receives the complete lambda arrays and its global state index; existing
temperature fields use the run temperature. Minimization has no added temperature
fields. `simulations_list.txt` lists window directory names in rest, coul, vdw order.
Window generation prints each window's global state index and bonded, Coulomb,
and VDW lambdas. The expected total charge is zero for a neutral input system;
this summary does not replace GROMACS preprocessing.

To customize the schedule, add campaign-level settings to the YAML used above:

```yaml
# Example custom schedule: 9 windows and 45 MDPs. Omit for the legacy 44 windows.
lambda_schedule:
  rest: [0.0, 0.5, 1.0]
  coul: [0.0, 0.5, 1.0]
  vdw: [0.0, 0.5, 1.0]

solvent_selection: "resname SOL NA CL B_NA B_CL"
```

Supply all three phase lists with the intended 0-to-1 endpoints. Values remain in
supplied order, including duplicates; counts and global indices follow those
lists. During rest, only bonded lambda changes; during coul, bonded lambda is 1
and Coulomb lambda changes; during vdw, bonded and Coulomb lambdas are 1 and VDW
lambda changes. Other lambda components are zero until their phase starts.
Repeated boundary states are retained. These optional settings do not change
`prepare_complex()`'s return value: the prepared complex directory.

The two steps can also be called independently on an existing alchemical complex:

```python
from vaasu.index import write_index_file
from vaasu.windows import create_windows

complex_path = "/data/DOR_10/abfe_van1_hrex_r1/complex"
index_path = write_index_file(complex_path)
windows = create_windows(complex_path, temperature=298.15)
```

`write_index_file()` returns one absolute path; `create_windows()` returns the
ordered list of absolute window paths. Existing index files, window directories,
or the simulation list are refused, including symlinks. Failures propagate and
leave partial output for inspection. The functions use supplied inputs directly;
there is no extra input-validation layer. TPR generation and simulation execution
are separate steps.

### ARCHER2 job scripts

`prepare_complex()` passes the paths returned by `create_windows()` directly to
`write_job_scripts()`. The returned order is used for every multidir stage. The
optional `job_name` argument defaults to the vanilla directory's parent name;
the YAML example above explicitly uses the ligand name.

VAASU writes the same five variants as the reference setup:

| Script | QoS | Walltime | Execution |
| --- | --- | --- | --- |
| `job_complex_archer_test.sh` | short | 00:20:00 | All stages |
| `job_complex_archer.sh` | standard | 24:00:00 | All stages |
| `job_complex_archer_contd.sh` | standard | 24:00:00 | Continue production |
| `job_complex_archer_lowpriority.sh` | lowpriority | 24:00:00 | All stages |
| `job_complex_archer_contd_lowpriority.sh` | lowpriority | 24:00:00 | Continue production |

The scripts retain the reference account `e280-biggin`, partition `standard`,
and GROMACS 2020.6 / PLUMED 2.8.2 installation paths under
`/work/e280/e280-Biggin/reub0138/gromacs_plumed_installation`. These cluster settings
are kept in the packaged `templates/jobs/archer.sh` template.

Resource requests use one MPI rank per window, up to two ranks per 128-core node,
and 64 CPUs per rank when there are at least two windows. Thus the default
44-window schedule uses 22 nodes; six windows use three nodes. An explicit total
task count handles an odd number of windows; a single window uses one node and
128 CPUs. See [Slurm task allocation options](https://slurm.schedmd.com/sbatch.html#OPT_ntasks-per-node).

Initial jobs preprocess each window on one rank, then run the stage together, in
`enmin`, `nvt`, `npt_b`, `npt_pr`, `prod` order. Every window uses
`complex_alchemical.top` and `index_alchemical.ndx`. Minimization starts from
`complex_alchemical.gro`; each later stage reads coordinates from its preceding
stage. Production retains `-replex 200 -hrex -dlb no` and uses
`plumed.dat` in the complex directory. Continuation jobs run only production with
`-cpi state.cpt`, using existing TPRs and checkpoints. Commands stop on failure.
An empty `plumed.dat` is created if absent; existing contents are preserved.

To generate jobs separately:

```python
from vaasu.jobs import write_job_scripts
from vaasu.windows import create_windows

windows = create_windows(complex_path, temperature=298.15)
job_paths = write_job_scripts(complex_path, windows=windows, job_name="DOR_10")
```

For a complex already prepared by `prepare_complex()`, the scripts already exist.
Copy the complete complex directory to ARCHER2 and submit from that directory:

```bash
sbatch job_complex_archer.sh
# After the production job has run, use its checkpoints to continue:
sbatch job_complex_archer_contd.sh
```

Job generation refuses existing script files or symlinks and preserves partial
output on failure. It generates the scripts locally without submitting jobs.

Existing `complex` or `boresch` paths raise `FileExistsError`.
Failures propagate and partial output is preserved. Restraint failures leave the
empty `complex/` directory and any partial `boresch/` output. Inspect and remove or
move partial directories explicitly before retrying. Source files are not changed.

After adding the two ITP arguments and YAML values, the Yuan example can be run
from any working directory after activating VAASU:

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
