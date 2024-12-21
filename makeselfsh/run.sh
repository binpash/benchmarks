#!/bin/bash

# Define repository structure and script paths
REPO_TOP=$(git rev-parse --show-toplevel)
testssl_dir="${REPO_TOP}/testssl"
scripts_dir="${testssl_dir}/source"
results_dir="${testssl_dir}/results"

# Define the target and output files
target="google.com" # Replace with the actual target or parameterize it
output_file="${results_dir}/testssl-output.txt"
error_file="${results_dir}/testssl-error.log"

# Create the results directory if it doesn't exist
mkdir -p "${results_dir}"

# Run the testssl.sh script and capture output
main_script="${scripts_dir}/testssl.sh"
"${main_script}" --json "${target}" > "${output_file}" 2> "${error_file}"

# Check for errors during execution
if [[ $? -ne 0 ]]; then
    echo "Error: testssl.sh script failed. Check '${error_file}' for details."
    exit 1
fi

echo "testssl.sh completed successfully. Output saved to '${output_file}'."

