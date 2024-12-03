#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional
import json
from subprocess import check_output, run
from collections import Counter
import os
from datetime import datetime

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

def get_parser():
    parser = argparse.ArgumentParser(
            prog='run_dynamic',
            description='runs the dynamic analysis')
    parser.add_argument('bench', type=str)
    parser.add_argument('--run-input', action=argparse.BooleanOptionalAction)
    parser.add_argument('--run-deps', action=argparse.BooleanOptionalAction)
    return parser

def get_environment(root):
    env = os.environ.copy()  
    dynamic_shell = root / 'infrastructure' / 'run_dynamic_shell.py'
    env['BENCHMARK_SHELL'] = str(dynamic_shell)
    env['BENCHMARK_EXPERIMENT_START'] = datetime.now().isoformat()
    return env

def run_analysis(root: Path, bench: Path, run_input: bool, run_deps: bool):
    env = get_environment(root)
    if run_deps:
        run([root / bench / 'deps.sh'], env=env)
    if run_input:
        run([root / bench / 'input.sh'], env=env)
    run([root / bench / 'run.sh'], env=env)

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    root = get_project_root()
    run_analysis(root, bench=args.bench, run_input=args.run_input, run_deps=args.run_deps)
