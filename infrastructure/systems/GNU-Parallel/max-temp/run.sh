#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

suffix=.full
if [[ "$@" == *"--small"* ]]; then
    suffix=.small
fi


export input_file="${input_dir}/temperatures$suffix.txt"
export statistics_dir="$results_dir/statistics$suffix"

mkdir -p "$statistics_dir"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="max-temp"
export BENCHMARK_INPUT_FILE="$(realpath "$input_file")"
export BENCHMARK_SCRIPT="$(realpath "${scripts_dir}/temp-analytics.sh")"

time $BENCHMARK_SHELL "${scripts_dir}/temp-analytics.sh"
