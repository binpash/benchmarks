#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
scripts_dir="${eval_dir}/scripts"
main_script="${scripts_dir}/git-workflow.sh"

export BENCHMARK_CATEGORY="git-workflow"
export BENCHMARK_SCRIPT="$main_script"
export BENCHMARK_INPUT_FILE="${eval_dir}/inputs/commits"

mkdir -p "${eval_dir}/outputs"

"$BENCHMARK_SHELL" "$main_script"
echo $?