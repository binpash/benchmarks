#!/usr/bin/env python3

import lzma
from collections import defaultdict
from pathlib import Path
import json
import math

from project_root import get_project_root

def correct_base(path):
    return Path(path).is_relative_to('/benchmarks')

def rebase(path):
    return Path(path).relative_to('/benchmarks')

def readings_dict(readings):
    return {r['log_current_time']: r for r in readings}

def is_shell(cmd):
#     non_shell = {'cat', 'tr', 'grep', 'sort', 'uniq', 'cut', 'awk', 'sed', 'rev', 'wc', 'convert', 'ffmpeg'}
#     is_shell_names = {'/bin/bash'}
#     assert cmd in non_shell | is_shell_names, f"unknown whether {cmd} is a shell"
#     return cmd in is_shell_names
    return cmd is not None and len(cmd) > 0 and 'bash' in cmd[0]

def sum_counters(pid, at_time, processes, children, should_include, which_counter):
    s = defaultdict(int)
    def recurse(pid):
        if at_time in processes[pid]:
            record = processes[pid][at_time]
            if should_include(record):
                for k, v in which_counter(record).items():
                    s[k] += v
        for c in children[pid]:
            recurse(c)
    recurse(pid)
    return s

def input_files(pid, at_time, processes, children):
    s = set()
    def recurse(pid):
        if at_time in processes[pid]:
            record = processes[pid][at_time]
            s.add(record['benchmark_input_file'])
            for file_path, _, _, mode, _ in record['full']['open_files']:
                known_modes = {'r', 'r+', 'w'}
                assert mode in known_modes, f"unknown mode {mode}"

                is_a_script = '.sh' in Path(file_path).suffixes
                if mode != 'w' and not is_a_script and correct_base(file_path):
                    s.add(file_path)
        for c in children[pid]:
            recurse(c)
    recurse(pid)
    s = s - {None}
    return s

def read_log_file(path):
    parents = defaultdict(lambda: None)
    children = defaultdict(set)
    processes = defaultdict(list)
    with lzma.open(path, 'r') as lines:
        for data in lines:
            data = json.loads(data)
            processes[data['pid']].append(data)
            pid = data['pid']
            parent = data['parent']
            children[parent].add(pid)
            parents[pid] = parent
    processes = defaultdict(lambda: None, {pid: readings_dict(rs) for pid, rs in processes.items()})
    return processes, parents, children

def print_statistics(pid, processes, parents, children):
    rs = processes[pid]

    max_uss = max(
        sum_counters(pid, log_time, processes, children, 
            lambda record: True,
            lambda record: record['pfullmem'],
        )
        ['uss'] 
        for log_time in rs
    )

    all_input_files = set()
    for log_time in rs: 
        all_input_files |= input_files(pid, log_time, processes, children)
    all_input_files = ";".join(str(rebase(p)) for p in all_input_files)

    max_reading = max(rs)
    tis = sum_counters(pid, max_reading, processes, children,
        lambda record: is_shell(record['cmdline']),
        lambda record: record['cpu_times'],
    )
    user = rs[max_reading]['cpu_times']['children_user']
    system = rs[max_reading]['cpu_times']['children_system']
    read_chars = rs[max_reading]['io_counters']['read_chars']
    write_chars = rs[max_reading]['io_counters']['write_chars']
    benchmark_script = rs[max_reading]['benchmark_script']
    benchmark_script = None if benchmark_script is None else rebase(benchmark_script)
    print(benchmark_script, user, system, max_uss, read_chars, write_chars, tis['user'], tis['system'], all_input_files, sep=',')

if __name__ == '__main__':
    process_logs = get_project_root() / 'infrastructure' / 'target' / 'process-logs'
    for path in process_logs.glob('*.jsonl.xz'):
        processes, parents, children = read_log_file(path)
        top_level = [pid for pid in processes if parents[parents[pid]] is None]
        for pid in top_level:
            print_statistics(pid, processes, parents, children)
