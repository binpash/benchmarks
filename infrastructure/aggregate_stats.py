#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from statistics import mean, stdev
from typing import Dict, List, Pattern

_PATTERNS: Dict[str, Pattern[str]] = {
    "Total CPU time [s]": re.compile(r"^Total CPU time:\s*([\d.]+)"),
    "Total Wall time [s]": re.compile(r"^Total Wall time:\s*([\d.]+)"),
    "Total IO bytes": re.compile(r"^Total IO bytes:\s*([\d.]+)"),
    "Max Memory Usage [bytes]": re.compile(r"^Max Memory Usage:\s*([\d.]+)"),
    "CPU time per input byte [s/byte]": re.compile(r"^CPU time per input byte:\s*([\d.]+)"),
    "Memory per input byte [bytes/byte]": re.compile(r"^Memory per input byte:\s*([\d.]+)"),
    "IO per input byte [bytes/byte]": re.compile(r"^IO per input byte:\s*([\d.]+)"),
    "Time in Shell [s]": re.compile(r"^Time in Shell:\s*([\d.]+)"),
    "Time in Commands [s]": re.compile(r"^Time in Commands:\s*([\d.]+)")
}

def _parse_file(path: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    with open(path, "r", encoding="utf‑8") as f:
        for line in f:
            for metric, pattern in _PATTERNS.items():
                m = pattern.match(line)
                if m:
                    out[metric] = float(m.group(1))
                    break
    return out


def _fmt(number: float | str) -> str:
    return f"{number:>12}"

def main(argv: List[str] | None = None) -> None:  # noqa: D401  (no period)
    parser = argparse.ArgumentParser(description="Aggregate benchmark statistics")
    parser.add_argument(
        "files",
        metavar="STAT_FILE",
        nargs="+",
        help="One or more *_bash_stats_run*.txt files to aggregate",
    )
    args = parser.parse_args(argv)

    aggregated: Dict[str, List[float]] = defaultdict(list)

    for path in args.files:
        try:
            stats = _parse_file(path)
        except OSError as exc:
            parser.error(f"Cannot read '{path}': {exc.strerror}")

        if not stats:
            print(f"Warning: no metrics found in '{path}', skipping", file=sys.stderr)
            continue

        for metric, value in stats.items():
            aggregated[metric].append(value)

    if not aggregated:
        parser.error("No metrics could be extracted from the supplied files")

    metric_width = max(len(m) for m in aggregated) + 2
    header = (
        f'{"Metric":<{metric_width}}'
        + _fmt("mean")
        + _fmt("min")
        + _fmt("max")
        + _fmt("stddev")
    )
    sep = "-" * len(header)

    print()
    print(f"Aggregated Statistics ({len(args.files)} files)")
    print(sep)
    print(header)
    print(sep)

    for metric in sorted(aggregated):
        values = aggregated[metric]
        m = mean(values)
        mn = min(values)
        mx = max(values)
        sd = stdev(values) if len(values) > 1 else 0.0
        print(f"{metric:<{metric_width}}" + _fmt(f"{m:.6f}") + _fmt(f"{mn:.6f}") + _fmt(f"{mx:.6f}") + _fmt(f"{sd:.6f}"))

    print(sep)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
