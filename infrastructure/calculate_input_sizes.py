#!/usr/bin/env python3
"""
generate_size_inputs.py  –  write size_inputs.jsonl for selected benchmarks.

Edit `BENCHMARKS` below to control which benchmark trees are scanned.
"""

import argparse
import json
import os
import pathlib
import subprocess
import sys
from typing import Iterable, List

BENCHMARKS: List[str] = ["aurpkg", "bio", "covid-mts", "dpt", "file-enc", "git-workflow", "llm", "log-analysis", "max-temp", "media-conv", "nlp", "oneliners", "prog-inf", "sklearn", "teraseq", "tuft-weather", "unix50"]
# riker is a special case

def git_root() -> pathlib.Path:
    try:
        root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
        return pathlib.Path(root)
    except subprocess.CalledProcessError:
        sys.exit("Not inside a Git repository.")


def du_size(path: str) -> int:
    try:
        out = subprocess.check_output(["du", "-b", path], text=True)
        return int(out.split()[0])
    except Exception:
        return os.stat(path).st_size


def walk_files(root: pathlib.Path):
    yield from (p for p in root.rglob("*") if p.is_file())


def emit_records(bench_root: pathlib.Path, out_file: pathlib.Path) -> None:
    with out_file.open("w", encoding="utf-8") as f:
        for bench in BENCHMARKS:
            inputs_dir = bench_root / bench / "inputs"
            if not inputs_dir.is_dir():
                sys.stderr.write(f"Skipped “{bench}” (no inputs/ directory)\n")
                continue

            for file in walk_files(inputs_dir):
                record = {
                    "size_bytes": du_size(str(file)),
                    "category": bench,
                    "path": file.relative_to(inputs_dir).as_posix(),
                }
                f.write(json.dumps(record, ensure_ascii=False) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate size_inputs.jsonl.")
    p.add_argument(
        "-r",
        "--root",
        default="benchmarks",
        help="Subdirectory (relative to repo root) that contains benchmarks "
             "(default: benchmarks/)",
    )
    p.add_argument(
        "-o",
        "--output",
        default="infrastructure/data/size_inputs.jsonl",
        help="Output JSONL file (default: size_inputs.jsonl)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    bench_root = git_root()

    emit_records(bench_root, pathlib.Path(args.output).resolve())


if __name__ == "__main__":
    main()
