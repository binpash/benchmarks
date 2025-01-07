#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/uniq-ips"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
cd "$eval_dir" # scripts/run.sh puts files in its current directory
$BENCHMARK_SHELL "$scripts_dir/run.sh" $@

