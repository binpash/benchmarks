#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"

export BENCHMARK_CATEGORY="riker"

for bench in "$scripts_dir"/*; do
    export BENCHMARK_SCRIPT="$(realpath "$bench/run.sh")"
    time "$bench/run.sh" "$@"
done

