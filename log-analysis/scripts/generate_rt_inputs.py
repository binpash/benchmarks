#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import random
from pathlib import Path
from typing import Iterable, Iterator

HEADER = "[RAY] pathID,timestamp,hop,value\n"
START_PATHID = 590_000              
BLOCK_LINES = 1_000_000            
NUM_FILES = 5
SEED = 42


def ray_line(pid: int, rng: random.Random) -> str:
    return (
        f"[RAY] {pid},{rng.randint(1000, 5000)},"
        f"{rng.randint(1, 20)},{rng.random():.4f}\n"
    )


def sequential_ids(start: int, count: int) -> Iterator[int]:
    for i in range(count):
        yield start + i


def random_ids(pool_start: int, pool_size: int, count: int,
               rng: random.Random) -> Iterator[int]:
    for _ in range(count):
        yield pool_start + rng.randint(0, pool_size - 1)


def write_info_file(
    path: Path,
    ids: Iterable[int],
    rng: random.Random,
    header: bool = False,
    block: int = BLOCK_LINES,
) -> None:
    with path.open("w") as fh:
        if header:
            fh.write(HEADER)

        buf: list[str] = []
        for pid in ids:
            buf.append(ray_line(pid, rng))
            if len(buf) >= block:
                fh.writelines(buf)
                buf.clear()
        if buf:
            fh.writelines(buf)


def generate_inputs(output_dir: str, num_lines: int, seed: int = SEED) -> None:
    out = Path(output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    rng = random.Random(seed)

    write_info_file(
        out / "1.INFO",
        sequential_ids(START_PATHID, num_lines),
        rng,
        header=True,
    )

    for idx in range(2, NUM_FILES + 1):
        write_info_file(
            out / f"{idx}.INFO",
            random_ids(START_PATHID, num_lines, num_lines, rng),
            rng,
            header=False,
        )

    print(
        f"Wrote {NUM_FILES} files × {num_lines:,} rows "
        f"→ {(NUM_FILES * num_lines):,} total lines under {out}"
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Create five *.INFO files for the ray-tracing benchmark"
    )
    ap.add_argument("output_dir", help="Directory in which to place the files")
    ap.add_argument("num_lines", type=int,
                    help="Number of data lines per file (N)")
    ap.add_argument("--seed", type=int, default=SEED,
                    help="RNG seed (default 42)")
    args = ap.parse_args()

    generate_inputs(args.output_dir, args.num_lines, args.seed)
