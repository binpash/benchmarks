#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/opt-parallel"
scripts_dir="${eval_dir}/scripts"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"

suffix=""

export input_dir="${eval_dir}/inputs/ChessData"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="opt-parallel"

main_script="${scripts_dir}/opt-parallel.sh"
BENCHMARK_SCRIPT="$(realpath "$main_script")"
export BENCHMARK_SCRIPT

"$BENCHMARK_SHELL" "$main_script" "$input_dir" > "${outputs_dir}/opt-parallel.out"
