"""Load ABFE campaign configuration from YAML."""

from pathlib import Path

import yaml


def load_config(path: str | Path) -> dict:
    """Read a campaign YAML file into ordinary dictionaries and lists."""
    with Path(path).open(encoding="utf-8") as config_file:
        return yaml.safe_load(config_file)
