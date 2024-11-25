#!/usr/bin/env python3

import shlex
from pathlib import Path
from typing import Optional
import json
from subprocess import check_output, run
from collections import Counter
import sys
from time import perf_counter

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

time_file = get_project_root() / 'infrastructure' / 'target' / 'runtime_log.csv'

command = sys.argv[1:]

start = perf_counter()
run(['bash', *command])
elapsed = perf_counter() - start

with time_file.open('a') as file:
    print(shlex.join(command), elapsed, sep=',', file=file)

