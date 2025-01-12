#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit-negate"
scripts_dir="${eval_dir}/scripts"
main_script="${scripts_dir}/vps-audit-negate.sh"
export BENCHMARK_CATEGORY=vps-audit-negate
export BENCHMARK_SCRIPT="$(realpath "$main_script")"
mkdir -p "${eval_dir}/outputs"
echo "Starting VPS audit..."
${BENCHMARK_SHELL} "${main_script}"
echo $?
