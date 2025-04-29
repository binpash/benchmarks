#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/calc/build.sh")"
export BENCHMARK_SCRIPT
(cd "$input_dir/calc/dev" && $BENCHMARK_SHELL "$scripts_dir/calc/build.sh")


