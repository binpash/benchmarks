#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="ci-cd"


for bench in "$eval_dir"/*; do
    if [ ! -d "$bench" ]; then
        continue
    fi
    export BENCHMARK_SCRIPT="$bench/execute.sh"
    $KOALA_SHELL "$bench/execute.sh" "$@"
done