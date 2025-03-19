#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
scripts_dir="${eval_dir}/scripts"
main_script_1="${scripts_dir}/vps-audit.sh"
main_script_2="${scripts_dir}/vps-audit-negate.sh"
export BENCHMARK_CATEGORY=vps-audit
export BENCHMARK_SCRIPT="$(realpath "$main_script_1")"
echo "Starting VPS audit..."
${BENCHMARK_SHELL} "${main_script_1}"
echo $?
export BENCHMARK_SCRIPT="$(realpath "$main_script_2")"
${BENCHMARK_SHELL} "${main_script_2}"
echo $?