#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
scripts_dir="${eval_dir}/scripts"
# Paths for the audit script and outputs
main_script="${scripts_dir}/vps-audit.sh"
error_file="${eval_dir}/outputs/vps-audit-error.log"
log_file="${eval_dir}/outputs/vps-audit-console.log"
output_file="${eval_dir}/outputs/vps-audit-output.txt"

# Create the outputs directory
mkdir -p "${eval_dir}/outputs"

# Run the audit script and capture real-time output in both log and output files
echo "Starting VPS audit..."
"${main_script}" 2>&1 | tee "${log_file}" > "${output_file}"
# Check for errors during execution
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "Error: VPS audit script failed. Check '${log_file}' for details."
    exit 1
fi

# Ensure any additional errors are captured
if [[ -s "${error_file}" ]]; then
    echo "Warnings or errors were logged during the audit. Check '${error_file}' for details."
fi

echo "VPS audit completed successfully. Output saved to '${output_file}' and console log to '${log_file}'."