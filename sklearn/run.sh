#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/sklearn"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
cd "$eval_dir" # scripts/run.sh references PWD
$BENCHMARK_SHELL "$scripts_dir/run.sh" $@

