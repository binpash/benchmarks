#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
scripts_dir="${eval_dir}/scripts"

# Execute the main audit script
main_script="${scripts_dir}/vps-audit.sh"

# Run the audit script and capture output
output_file="${eval_dir}/results/vps-audit-output.txt"
error_file="${eval_dir}/results/vps-audit-error.log"

mkdir -p "${eval_dir}/results"

"${main_script}" --output "${output_file}" > "${output_file}" 2> "${error_file}"

# Check for errors during execution
if [[ $? -ne 0 ]]; then
    echo "Error: VPS audit script failed. Check '${error_file}' for details."
    exit 1
fi

echo "VPS audit completed successfully. Output saved to '${output_file}'."