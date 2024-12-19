#!/usr/env/bin python3

from collections import defaultdict
import json
import math

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

path = get_project_root() / 'infrastructure' / 'target' / 'process_data_nlp.jsonl'

def sum_counters(p, which_counter, s = defaultdict(lambda: -math.inf)):
    for key, value in p[which_counter].items():
        s[key] += value
    for c in p['children']:
        sum_counters(c, which_counter, s)
    return s

with path.open('r') as data:

    processes = defaultdict(list)
    for data in data:
        data = json.loads(data)
        processes[data['pid']].append(data)

    processes = {pid: {'children': [], 'parent': None, 'readings': readings} for pid, readings in processes.items()}
    for pid, p in processes.items():
        q = p['readings'][0]
        parent = q['parent']
        if parent in processes:
            p['parent'] = processes[parent]
    for pid, p in processes.items():
        parent = p['parent']
        if parent is not None:
            parent['children'].append(p)
    processes = {k: v for k, v in processes.items() if v['parent'] is None}

    for p in processes.values():
        defaultdict(lambda: -math.inf)
        user = -math.inf
        system = -math.inf
        for r in p['readings']:
            user = max(user, r['cpu_times']['children_user'])
            system = max(system, r['cpu_times']['children_system'])

        print(user, system, p['readings'][0]['benchmark_script'])

