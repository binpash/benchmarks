#!/usr/bin/env python3

from subprocess import run, CalledProcessError
from pathlib import Path

def get_project_root():
    result = run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(f'could not find project root: `{result.stderr}`')
    return Path(result.stdout.removesuffix('\n')) # git only emits one trailing newline in the path
