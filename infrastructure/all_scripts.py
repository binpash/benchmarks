#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json

from project_root import get_project_root

def get_all_scripts(
    scripts_file: Path = get_project_root() / 'infrastructure/data/script-globs.json'
) -> list[Path]:
    scripts = scripts_file.read_text()
    benchmark_data: dict[str, dict[str, any]] = json.loads(scripts)
    return {
        benchmark_name: [
            script
            for script_glob in benchmark_data['scripts']
            for script in get_project_root().glob(script_glob)
        ]
        for benchmark_name, benchmark_data in benchmark_data.items()
    }

if __name__ == "__main__":
    for bench in get_all_scripts().keys():
        print(bench)
    print(get_all_scripts())
