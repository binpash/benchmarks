#!/usr/bin/env python3

import lzma
import tempfile
import argparse
from pathlib import Path
from typing import Optional
import json
from subprocess import check_output, run
from collections import Counter
import os
from datetime import datetime, timezone

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root


def get_parser():
    parser = argparse.ArgumentParser(
        prog="run_dynamic",
        description="runs the dynamic analysis",
    )
    parser.add_argument("bench", type=str)
    parser.add_argument(
        "forward", nargs=argparse.REMAINDER, help="pass these arguments to execute.sh"
    )
    return parser


def get_environment(
    root: Path, start_time: str, bench: str, data_log: str, mortem_log: str
):
    env = os.environ.copy()
    dynamic_shell = root / "infrastructure" / "run_dynamic_shell.py"
    io_shell = root / "infrastructure" / "io_shell.py"
    env["BENCHMARK_IO_SHELL"] = str(io_shell)
    env["KOALA_SHELL"] = str(dynamic_shell)
    env["BENCHMARK_EXPERIMENT_START"] = start_time
    env["BENCHMARK_PROCESS_LOG"] = data_log
    env["BENCHMARK_MORTEM_LOG"] = mortem_log
    return env


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    root = get_project_root()
    raw_bench = args.bench
    bench = Path(raw_bench.rstrip("/")).name
    start_time = datetime.now(timezone.utc).isoformat()
    bench_id = f"{start_time}" f"-" f'{str(bench).replace("/", "%")}'

    bench_dir = root / bench
    execute = bench_dir / "execute.sh"
    output_file = bench_dir / f"{bench}.out"
    error_file = bench_dir / f"{bench}.err"

    with tempfile.NamedTemporaryFile() as data_log:
        data_log = data_log.name
        mortem_log = (
            root / "infrastructure" / "target" / "process-logs" / f"{bench_id}.mortem"
        )
        env = get_environment(
            root=root,
            start_time=start_time,
            bench=bench,
            data_log=data_log,
            mortem_log=mortem_log,
        )
        # write to an uncompressed file because it is faster
        # run([root / bench / 'execute.sh', *args.forward], env=env, cwd=root / bench)
        output_file = root / bench / f"{bench}.out"
        error_file = root / bench / f"{bench}.err"

        with open(output_file, "w") as out, open(error_file, "w") as err:
            run(
                [root / bench / "execute.sh", *args.forward],
                env=env,
                cwd=root / bench,
                stdout=out,
                stderr=err,
            )

        compressed_data_log = (
            root / "infrastructure" / "target" / "process-logs" / f"{bench_id}.jsonl.xz"
        )
        with compressed_data_log.open("w") as stdout:
            run(["xz", "-7e", "-T0", "-c", data_log], stdout=stdout)
