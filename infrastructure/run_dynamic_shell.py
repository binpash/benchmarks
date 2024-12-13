#!/usr/bin/env python3

import os
from datetime import datetime
from itertools import chain
import time
import psutil
import signal
from pathlib import Path
from typing import Optional
import json
from subprocess import run
import sys
import asyncio

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

def data_json(p: psutil.Process) -> str: 
    p = p.as_dict()
    times = p['cpu_times']
    mem = p['memory_full_info']
    io_counters = p['io_counters']
    return json.dumps({
        'pid': p['pid'],
        'benchmark_category': p['environ'].get('BENCHMARK_CATEGORY'),
        'benchmark_script': p['environ'].get('BENCHMARK_SCRIPT'),
        'benchmark_input_file': p['environ'].get('BENCHMARK_INPUT_FILE'),
        'benchmark_experiment_start': p['environ'].get('BENCHMARK_EXPERIMENT_START'),
        'log_current_time': datetime.now().isoformat(),
        'cwd': p['cwd'],
        'cmdline': p['cmdline'],
        'create_time': p['create_time'],
        'cpu_times': {
            'user': times.user,
            'system': times.system,
            'children_user': times.children_user,
            'children_system': times.children_system,
            'iowait': times.iowait,
        },
        'pfullmem': {
            'rss': mem.rss, 'vms': mem.vms, 'shared': mem.shared, 'text': mem.text, 'lib': mem.lib, 'data': mem.data, 'dirty': mem.dirty, 'uss': mem.uss, 'pss': mem.pss, 'swap': mem.swap,
        },
        'io_counters': {
            'read_count': io_counters.read_count, 'write_count': io_counters.write_count, 'read_bytes': io_counters.read_bytes, 'write_bytes': io_counters.write_bytes, 'read_chars': io_counters.read_chars, 'write_chars': io_counters.write_chars,
        },
        'num_fds': p['num_fds'],
        'full': p, # this does not provide field names for cpu_times and io_counters, etc.
    })

def write_process_data(parent: int, data_log):
    parent = psutil.Process(parent)
    for p in chain(parent.children(recursive=True), [parent]):
        try:
            print(data_json(p), file=data_log)
        except psutil.NoSuchProcess:
            pass

async def collect_process_data(parent: int, stdout):
    try:
        write_process_data(parent, stdout)
        while True:
            await asyncio.sleep(0.2)
            write_process_data(parent, stdout)
    except Exception as e:
        print(e, type(e))

async def run_and_collect(program, data_log: Path):
    start_time = time.perf_counter()
    process = await asyncio.create_subprocess_exec(*program)
    pid = process.pid
    with data_log.open('a') as data_log:
        process_data = asyncio.create_task(collect_process_data(pid, data_log))
        await process.wait()
        print('pid', pid, time.perf_counter() - start_time)
        process_data.cancel()

async def main():
    program = sys.argv[1:]
    category = os.environ.get('BENCHMARK_CATEGORY')
    data_log = get_project_root() / 'infrastructure' / 'target' / f'process_data_{category}.jsonl'
    await run_and_collect(program=program, data_log=data_log)

asyncio.run(main())
