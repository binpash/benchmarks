#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
scripts_dir="${eval_dir}/scripts"
main_script="${scripts_dir}/vps-audit.sh"
log_file="${eval_dir}/outputs/vps-audit-console.log"

mkdir -p "${eval_dir}/outputs"

# run the audit script and capture real-time output
echo "Starting VPS audit..."
"${main_script}" 2>&1 | tee "${log_file}"
# Check for errors during execution
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "Error: VPS audit script failed. Check '${log_file}' for details."
    exit 1
fi

echo "VPS audit completed successfully. Output saved to '${output_file}' and console log to '${log_file}'."