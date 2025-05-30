#!/bin/bash

# For sysctl
export PATH="$PATH:/sbin:/usr/sbin"

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
scripts_dir="${eval_dir}/scripts"
main_script_1="${scripts_dir}/vps-audit.sh"
main_script_2="${scripts_dir}/vps-audit-negate.sh"
export BENCHMARK_CATEGORY=vps-audit
BENCHMARK_SCRIPT="$(realpath "$main_script_1")"
export BENCHMARK_SCRIPT
echo "Starting VPS audit..."
${KOALA_SHELL} "${main_script_1}"
echo $?
BENCHMARK_SCRIPT="$(realpath "$main_script_2")"
export BENCHMARK_SCRIPT
${KOALA_SHELL} "${main_script_2}"
echo $?
