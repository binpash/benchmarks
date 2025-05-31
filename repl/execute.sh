#!/bin/bash

# For sysctl
export PATH="$PATH:/sbin:/usr/sbin"

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/repl"
scripts_dir="${eval_dir}/scripts"
main_script_1="${scripts_dir}/vps-audit.sh"
main_script_2="${scripts_dir}/vps-audit-negate.sh"

export BENCHMARK_CATEGORY=repl

BENCHMARK_SCRIPT="$(realpath "$main_script_1")"
export BENCHMARK_SCRIPT
echo "Starting VPS audit..."
${KOALA_SHELL} "${main_script_1}"
echo $?

BENCHMARK_SCRIPT="$(realpath "$main_script_2")"
export BENCHMARK_SCRIPT
${KOALA_SHELL} "${main_script_2}"
echo $?

git_script="${scripts_dir}/git-workflow.sh"

export BENCHMARK_SCRIPT="$git_script"
export BENCHMARK_INPUT_FILE="${eval_dir}/inputs"

mkdir -p "${eval_dir}/outputs"

NUM_COMMITS=21

for arg in "$@"; do
    case "$arg" in
        --min) NUM_COMMITS=2 ;;
        --small) NUM_COMMITS=6 ;;
    esac
done
$KOALA_SHELL "$git_script" "$NUM_COMMITS" "$@"
echo $?
