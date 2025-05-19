#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/llm"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="llm"


for bench in "$scripts_dir"/*; do
    export BENCHMARK_SCRIPT="$bench/scripts/execute.sh"
    $KOALA_SHELL "$bench/execute.sh" "$@"
done