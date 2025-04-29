#!/usr/bin/env python3
"""
aggregate_stats.py  FILE1 FILE2 ...

Reads one or more benchmark_stats_run*.txt files and prints averages,
minima and maxima for each numeric metric found.
"""

import re, sys, statistics, pathlib, textwrap

num_re = re.compile(r"^(.*?):\s+([0-9.]+)")
wanted = {
    "Total CPU time",
    "Total Wall time",
    "Total IO bytes",
    "Max Memory Usage",
    "CPU time per input byte",
    "Memory per input byte",
    "IO per input byte",
    "Time in Shell",
    "Time in Commands",
}

metrics = {key: [] for key in wanted}

def parse_file(path: pathlib.Path):
    with path.open() as fh:
        for line in fh:
            m = num_re.match(line.strip())
            if m and m[1] in wanted:
                metrics[m[1]].append(float(m[2]))

for p in map(pathlib.Path, sys.argv[1:]):
    if not p.exists():
        sys.exit(f"File not found: {p}")
    parse_file(p)

print("Aggregated Benchmark Statistics")
print("=" * 42)
print(f"Runs analysed: {len(sys.argv) - 1}\n")

for k, vals in metrics.items():
    if not vals:     
        continue
    avg = statistics.mean(vals)
    mn, mx = min(vals), max(vals)
    print(textwrap.dedent(f"""\
        {k}:
            mean = {avg:,.6f}
            min  = {mn:,.6f}
            max  = {mx:,.6f}
    """))
