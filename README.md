# VAASU

VAASU's An ABFEP Setup Utility: a Python library with its own uv environment.
ABFE tools are not implemented yet. Requires uv and Python 3.12.

## Setup and activation

```bash
cd /biggin/b211/reub0138/Projects/VAASU
uv sync --locked
source .venv/bin/activate
```

This creates `.venv` and installs VAASU in editable mode with Ruff.
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

## Development

Run from the repository; replace `<package>` with the dependency name.

```bash
uv add <package>          # Runtime dependency
uv add --dev <package>    # Development dependency
uv remove <package>
uv run ruff check .
uv run ruff format .
uv build                 # Wheel and source distribution in dist/
```

Commit `pyproject.toml` and `uv.lock` together after dependency changes.
