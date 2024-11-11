#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json
from subprocess import check_output

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes

for benchmark_name, scripts in get_all_scripts().items():
    total = 0
    for script in scripts:
        cloc = check_output(['cloc', '--json', script], text=True)
        cloc = json.loads(cloc)['SUM']['code']
        total += cloc
#         print(benchmark_name, script, cloc)
    print('cloc', benchmark_name, script, total)

for benchmark_name, scripts in get_all_scripts().items():
    for script in scripts:
        asts = parse_shell_script(script)
        print('node count', benchmark_name, script, count_nodes(asts))
