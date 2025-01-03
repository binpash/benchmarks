#!/usr/bin/env python3

import os
from datetime import datetime, timezone
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

def data_json(p: psutil.Process, log_current_time: str, benchmark_experiment_start: str, benchmark_script: Path) -> str: 
    parent_pid = None
    try:
        parent_pid = p.parent().pid
    except AttributeError:
        pass
    p = p.as_dict()
    times = p['cpu_times']
    mem = p['memory_full_info']
    io_counters = p['io_counters']
    return json.dumps({
        'pid': p['pid'],
        'parent': parent_pid,
        'benchmark_category': p['environ'].get('BENCHMARK_CATEGORY'),
        'benchmark_script': str(benchmark_script),
        'benchmark_input_file': p['environ'].get('BENCHMARK_INPUT_FILE'),
        'benchmark_experiment_start': benchmark_experiment_start,
        'log_current_time': log_current_time,
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

def write_process_data(parent: int, data_log, benchmark_experiment_start, benchmark_script: Path):
    log_current_time = datetime.now(timezone.utc).isoformat()
    parent = psutil.Process(parent)
    for p in chain(parent.children(recursive=True), [parent]):
        try:
            print(data_json(p, log_current_time, benchmark_experiment_start, benchmark_script), file=data_log)
        except psutil.NoSuchProcess:
            pass

async def collect_process_data(parent: int, data_log, benchmark_experiment_start, benchmark_script):
    try:
        write_process_data(parent, data_log, benchmark_experiment_start, benchmark_script)
        while True:
            await asyncio.sleep(0.05)
            write_process_data(parent, data_log, benchmark_experiment_start, benchmark_script)
    except Exception as e:
        print(e, type(e), file=sys.stderr)

async def run_and_collect(
    program, data_log: Path, mortem_log: Path, benchmark_experiment_start: Path, root: Path, benchmark_script: Path, 
):
    start_time = time.perf_counter()
    process = await asyncio.create_subprocess_exec(*program)
    pid = process.pid
    with data_log.open('a') as stdout:
        process_data = asyncio.create_task(collect_process_data(pid, stdout, benchmark_experiment_start, benchmark_script))
        await process.wait()
        end_time = time.perf_counter() 
        process_data.cancel()
    with mortem_log.open('a') as mortem_log:
        print(benchmark_script.relative_to(root), benchmark_experiment_start, pid, end_time - start_time, sep=',', file=mortem_log)

async def main():
    program = sys.argv[1:]
    category = os.environ.get('BENCHMARK_CATEGORY')
    data_log = Path(os.environ.get('BENCHMARK_PROCESS_LOG'))
    mortem_log = Path(os.environ.get('BENCHMARK_MORTEM_LOG'))
    benchmark_experiment_start = os.environ.get('BENCHMARK_EXPERIMENT_START')
    benchmark_script = Path(os.environ.get('BENCHMARK_SCRIPT'))
    root = Path(get_project_root())
    await run_and_collect(
        program=program, 
        data_log=data_log, 
        mortem_log=mortem_log,
        benchmark_experiment_start=benchmark_experiment_start,
        root=root,
        benchmark_script=benchmark_script,
    )

asyncio.run(main())
