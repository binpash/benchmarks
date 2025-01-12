#!/usr/bin/env python3

from time import perf_counter
from sys import argv

duration = float(argv[1])
start = perf_counter()
while perf_counter() - start < duration:
    for _ in range(int(1e6)):
        pass
