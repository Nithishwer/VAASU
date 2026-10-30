# VAASU

VAASU stands for VAASU's An ABFEP Setup Utility.

VAASU is a Python library with its own uv-managed environment, independent of
abflow. The package foundation is in place; ABFE tools and scientific dependencies
will be added as they are implemented. The old ABFE setup script can use this
environment once its required VAASU tools are available.

## Set up the environment

Use uv and Python 3.12. From the repository directory, run:

```bash
cd /biggin/b211/reub0138/Projects/VAASU
uv sync --locked
```

This creates `.venv`, installs the dependencies recorded in `uv.lock`, and installs
VAASU in editable mode. Changes under `src/vaasu` are then available without
reinstalling the package. Ruff is included as a development tool.

Run `uv sync --locked` again after pulling dependency changes. The environment is
local to this checkout and should not be committed to Git.

Python 3.12 is the initial development baseline. The modified restraint generator
will be validated when it is introduced; Python 3.11 will be considered if a
required dependency demonstrates an incompatibility.

## Activate and deactivate the environment

In a Bash terminal, activate the environment:

```bash
source /biggin/b211/reub0138/Projects/VAASU/.venv/bin/activate
```

Activation makes the environment's Python available directly, including when you
change to another directory:

```bash
python
python /absolute/path/to/script.py
```

Exit the interactive Python session before running further shell commands. To
check which interpreter is selected, run:

```bash
python -c "import sys; print(sys.executable)"
```

It should print `/biggin/b211/reub0138/Projects/VAASU/.venv/bin/python`.

When finished, leave the environment:

```bash
deactivate
```

## Run without activating

From the repository directory, uv selects the project environment automatically:

```bash
uv run python
uv run python /absolute/path/to/script.py
```

From another directory, select the project explicitly:

```bash
uv run --project /biggin/b211/reub0138/Projects/VAASU python /absolute/path/to/script.py
```

Replace `/absolute/path/to/script.py` with your script's path. These commands use
VAASU's environment without changing your shell's active environment. Scripts can
use `import vaasu` without modifying `sys.path` or `PYTHONPATH`.

## Manage dependencies

Run these commands from the repository directory, replacing `<package>` with the
package name:

```bash
uv add <package>
uv add --dev <package>
uv remove <package>
```

Use `uv add` for a dependency needed by the library and `uv add --dev` for a
development tool. To remove a development dependency, use
`uv remove --dev <package>`.

These commands update the dependency declarations, lockfile, and environment.
Commit `pyproject.toml` and `uv.lock` together after dependency changes. Add
scientific dependencies only when an implemented tool needs them.

## Check, format, and build

From the repository directory, check the code and apply Python formatting:

```bash
uv run ruff check .
uv run ruff format .
```

To check formatting without changing files:

```bash
uv run ruff format --check .
```

Build an installable wheel and source distribution:

```bash
uv build
```

The artifacts are written to `dist/`, which is excluded from Git.
