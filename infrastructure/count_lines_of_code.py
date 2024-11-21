#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json
from subprocess import Popen, PIPE

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

root = get_project_root()
for benchmark_name, scripts in get_all_scripts().items():
    processes = []
    for script in scripts:
        process = Popen(['cloc', '--json', script], stdout=PIPE)
        script = script.relative_to(root)
        processes.append((script, process))
    for script, process in processes:
        stdout, _stderr = process.communicate()
        stdout = stdout.decode()
        cloc = json.loads(stdout)
        cloc = cloc['SUM']['code']
        print(script, cloc, sep=',')
