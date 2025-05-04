#!/bin/bash

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
scripts_dir="${eval_dir}/scripts"
main_script="${scripts_dir}/git-workflow.sh"

export BENCHMARK_CATEGORY="git-workflow"
export BENCHMARK_SCRIPT="$main_script"
export BENCHMARK_INPUT_FILE="${eval_dir}/inputs"

mkdir -p "${eval_dir}/outputs"

NUM_COMMITS=21

for arg in "$@"; do
    case "$arg" in
        --min) NUM_COMMITS=2 ;;
        --small) NUM_COMMITS=6 ;;
    esac
done
$KOALA_SHELL "$main_script" "$NUM_COMMITS" "$@"
echo $?