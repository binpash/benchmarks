#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json
from subprocess import check_output

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

root = get_project_root()
for benchmark_name, scripts in get_all_scripts().items():
    for script in scripts:
        print(script.relative_to(root), benchmark_name, sep=',')
