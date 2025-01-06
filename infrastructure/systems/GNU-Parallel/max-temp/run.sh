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

$BENCHMARK_SHELL ${scripts_dir}/temp-analytics.sh
