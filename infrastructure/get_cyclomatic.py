#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
from typing import Optional
import json
from subprocess import check_output

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

root = get_project_root()
shellmetrics = root / 'infrastructure' / 'target' / 'shellmetrics.sh'
if not shellmetrics.is_file():
    raise f'shellmetrics.sh not in `{shellmetrics}`, please put the script at that path'

all_scripts = [s for scripts in get_all_scripts().values() for s in scripts]

output = check_output([shellmetrics, '--csv', '--shell', 'bash', '--no-color', *all_scripts], text=True)
datas = defaultdict(list)
for line in output.splitlines()[1:]:
    file, _func, _lineno, _lloc, ccn, _lines, _comment, _blank = line.split(',')
    file = json.loads(file)
    file = Path(file).relative_to(root)
    datas[file].append((ccn,))

for file, datas in datas.items():
    ccn = sum(float(ccn) for ccn, in datas)
    print(file, ccn, sep=',')
