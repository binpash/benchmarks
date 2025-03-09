#!/usr/bin/env python3

from pathlib import Path
from subprocess import Popen
from time import perf_counter, sleep
from os import getpid, environ
from sys import argv

# need zombie io from procfs 
# why? because it will include exactly the io the process rchar and wchar

# need zombie stat from procfs 
# why? because it will include user time and system time (exactly time in shell)

# need stat from procfs before and after process exits,
# why? because it will include cutime and cstime (subtract to get total user and system time)

stat_path = Path('/proc', str(getpid()), 'stat')

before_stat = stat_path.read_bytes()
before_seconds = perf_counter()

# must execute actual benchmark script here because we use its procfs utime and stime for the time in shell
with Popen(argv[1:]) as proc:
    pid = str(proc.pid)
    io_zombie_path = Path('/proc', pid, 'io')
    stat_zombie_path = Path('/proc', pid, 'stat')
    while True:
        # the last time these are set, they are collected from a zombie (probably?)
        # we need admin privilages to collect the data from a zombie
        io_zombie = io_zombie_path.read_bytes()
        stat_zombie = stat_zombie_path.read_bytes()
        if proc.poll() is not None:
            break
        sleep(1e-1) # assume process exits during this call
after_seconds = perf_counter() # may miss up to 0.1 seconds
after_stat = stat_path.read_bytes()

elapsed = after_seconds - before_seconds
Path(environ['BENCHMARK_ELAPSED_SECS']).write_text(str(elapsed))
Path(environ['BENCHMARK_IO_ZOMBIE']).write_bytes(io_zombie)
Path(environ['BENCHMARK_STAT_ZOMBIE']).write_bytes(stat_zombie)
Path(environ['BENCHMARK_STAT_BEFORE']).write_bytes(before_stat)
Path(environ['BENCHMARK_STAT_AFTER']).write_bytes(after_stat)