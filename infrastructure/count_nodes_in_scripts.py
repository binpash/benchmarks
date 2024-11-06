#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes

scripts = get_all_scripts()

for script in scripts:
    asts = parse_shell_script(script)
    print(str(script), count_nodes(asts))
