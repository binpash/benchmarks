#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/sklearn"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="sklearn"
 
cd "$eval_dir" || exit 1 # scripts/execute.sh references PWD
 
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/execute.sh")"
export BENCHMARK_SCRIPT
$KOALA_SHELL "$scripts_dir/execute.sh" "$@"

