#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/uniq-ips"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="uniq-ips"

cd "$eval_dir" # scripts/run.sh puts files in its current directory

export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/run.sh")"
$BENCHMARK_SHELL "$scripts_dir/run.sh" $@

