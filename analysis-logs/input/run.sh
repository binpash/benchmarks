#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-logs"
results_dir="${eval_dir}/results"
inputs_dir="${eval_dir}/input"

shell="/bin/bash"

mkdir -p $results_dir

export INPUT=${inputs_dir}/access.log
script="${eval_dir}/nginx.sh"

echo "Executing $(basename "$script")"
$shell "$script" > "$results_dir/$(basename "$script").out"
