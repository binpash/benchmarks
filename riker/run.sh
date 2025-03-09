#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="riker"

for bench in "$scripts_dir"/*; do
    "$BENCHMARK_SHELL" "$bench/run.sh" "$@"
done