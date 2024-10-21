#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-logs"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"
mkdir -p $results_dir

export INPUT=${input_dir}/access.log

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

log_dir="$results_dir/results$suffix"
mkdir -p $log_dir
$scripts_dir/nginx.sh > $log_dir/out
