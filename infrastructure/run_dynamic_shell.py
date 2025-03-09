#!/usr/bin/env python3

from tempfile import TemporaryDirectory
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

def get_data_json(
    p: psutil.Process, 
    log_current_time: str, 
    benchmark_experiment_start: str, 
) -> str: 
    parent_pid = p.ppid()
    mem = p.memory_full_info() # slow, we need uss
    return dict(
        pid=p.pid,
        parent=parent_pid,
        benchmark_experiment_start=benchmark_experiment_start, # identifies top level script run
        log_current_time= log_current_time, # identifies log
        cwd= p.cwd(),
        times=list(p.cpu_times()),
        cmdline= p.cmdline(),
        create_time= p.create_time(),
        uss= mem.uss,
        num_fds= p.num_fds(),
        open_files= p.open_files(),
    )

def write_process_data(parent: int, data_log, benchmark_experiment_start):
    parent = psutil.Process(parent)
    log_current_time = datetime.now(timezone.utc).isoformat()
    procs = [parent, *parent.children(recursive=True)]
    logs = []
    for p in procs:
        try: # exception most likely that we can't read procfs because process doesn't exist. collect what we can.
            data = get_data_json(p, log_current_time, benchmark_experiment_start)
            logs.append(data)
        except psutil.NoSuchProcess: # most likely
            pass
        except Exception as e:
            print(e, type(e), file=sys.stderr)
            continue
    data_log.append(logs)

async def collect_process_data(parent, data_log, benchmark_experiment_start: str):
    write_process_data(parent, data_log, benchmark_experiment_start)
    while True:
        await asyncio.sleep(1e-2)
        try: 
            write_process_data(parent, data_log, benchmark_experiment_start)
        except psutil.NoSuchProcess: # most likely
            pass
        except Exception as e:
            print(e, type(e), file=sys.stderr)

io_shell_data_items = ['io_zombie', 'stat_zombie', 'elapsed_secs', 'stat_before', 'stat_after']
def get_environment(io_shell_data_dir: Path):
    env = os.environ.copy()
    for item in io_shell_data_items:
        env[f'BENCHMARK_{item.upper()}'] = str(io_shell_data_dir / item)
    return env

def get_io_shell_data(io_shell_data_dir: Path):
    return {item: (io_shell_data_dir / item).read_text() for item in io_shell_data_items}

async def main():
    root = Path(get_project_root())
    program = sys.argv[1:]
    io_shell = Path(os.environ['BENCHMARK_IO_SHELL'])
    category = os.environ['BENCHMARK_CATEGORY']
    # BENCHMARK_INPUT_FILE is optional: it is not trivial to get it for all benchmarks
    input_file = os.environ.get('BENCHMARK_INPUT_FILE')
    script = os.environ['BENCHMARK_SCRIPT']
    benchmark_experiment_start = os.environ['BENCHMARK_EXPERIMENT_START']
    data_log = Path(os.environ['BENCHMARK_PROCESS_LOG'])
    mortem_log = Path(os.environ['BENCHMARK_MORTEM_LOG'])
    benchmark_script = Path(os.environ['BENCHMARK_SCRIPT'])
    data_log_list = []
    with TemporaryDirectory() as io_shell_data_dir:
        io_shell_data_dir = Path(io_shell_data_dir)
        env = get_environment(io_shell_data_dir)
        process = await asyncio.create_subprocess_exec(io_shell, *program, env=env)
        pid = process.pid
        process_data_task = asyncio.create_task(collect_process_data(
            pid, 
            data_log_list, 
            benchmark_experiment_start, 
        ))
        await process.wait()
        process_data_task.cancel()
        try:  
            await process_data_task
        except asyncio.exceptions.CancelledError:
            pass
        log_entry = get_io_shell_data(io_shell_data_dir)
    log_entry |= dict(
        script=script,
        pid=pid,
        benchmark_experiment_start=benchmark_experiment_start,
        category=category,
        input_file=input_file,
        sc_clk_tck=os.sysconf('SC_CLK_TCK'),
    )
    with mortem_log.open('a') as mortem_log:
        json.dump(log_entry, mortem_log)
        print(file=mortem_log)
    with data_log.open('a') as data_log:
        for log_entry in data_log_list:
            json.dump(log_entry, data_log)
            print(file=data_log)

asyncio.run(main())