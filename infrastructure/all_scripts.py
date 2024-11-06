#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json

from project_root import get_project_root

def get_all_scripts(
    scripts_file: Path = get_project_root() / 'infrastructure/data/script-globs.json'
) -> list[Path]:
    scripts = scripts_file.read_text()
    script_globs: list[str] = json.loads(scripts)
    return sorted(
        script
        for script_glob in script_globs
        for script in get_project_root().glob(script_glob)
    )
