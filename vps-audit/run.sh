#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
scripts_dir="${eval_dir}/scripts"
script_1="${scripts_dir}/vps-audit.sh"
script_2="${scripts_dir}/vps-audit-negate.sh"
mkdir -p "${eval_dir}/outputs"
echo "Starting VPS audit..."
${BENCHMARK_SHELL} "${script_1}" > "${eval_dir}/outputs/vps-audit.log"
echo $?
echo "Starting VPS audit negate..."
${BENCHMARK_SHELL} "${script_2}" > "${eval_dir}/outputs/vps-audit-negate.log"
echo $?
