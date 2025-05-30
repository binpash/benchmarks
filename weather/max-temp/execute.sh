#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/inputs"

export LC_ALL=C

suffix=.full
for arg in "$@"; do
    if [ "$arg" = "--small" ]; then
        suffix=".small"
        break
    elif [ "$arg" = "--min" ]; then
        suffix=".min"
        break
    fi
done


export input_file="${input_dir}/temperatures$suffix.txt"
export statistics_dir="$outputs_dir/statistics$suffix"

mkdir -p "$statistics_dir"

KOALA_SHELL=${KOALA_SHELL:-bash}

BENCHMARK_CATEGORY="max-temp"
export BENCHMARK_CATEGORY

BENCHMARK_INPUT_FILE="$(realpath "$input_file")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "${scripts_dir}/temp-analytics.sh")"
export BENCHMARK_SCRIPT

$KOALA_SHELL "$scripts_dir/temp-analytics.sh"
