#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

shell="/bin/bash"

mkdir -p $results_dir

export input_file=${input_dir}/temperatures2015.txt
export results_dir
script="${scripts_dir}/temp-analytics.sh"

echo "Executing $(basename "$script")"
$shell "$script" > "$results_dir/$(basename "$script").out"
