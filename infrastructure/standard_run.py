#!/usr/bin/env python3

from pathlib import Path
from subprocess import run

from project_root import get_project_root

root = get_project_root()
script_dir = Path(__file__).resolve().parent
standard_dir = script_dir / 'standards'

for s in sorted(standard_dir.iterdir()):
    print('running standard', s.relative_to(root))
    run([s / 'deps.sh'])
    run([s / 'input.sh'])
    run([script_dir / 'run_dynamic.py', s.relative_to(root)])
    run([s / 'cleanup.sh'])


