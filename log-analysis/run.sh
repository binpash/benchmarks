#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/inputs"
scripts_dir="${eval_dir}/scripts"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"

nginx_input=$input_dir/nginx-logs
pcaps_input=$input_dir/pcaps
suffix=".full"
for arg in "$@"; do
    if [ "$arg" = "--small" ]; then
        # TODO: vary input size
        nginx_input=$input_dir/nginx-logs
        pcaps_input=$input_dir/pcaps
        suffix=".small"
    fi
    if [ "$arg" = "--min" ]; then
        nginx_input=$input_dir/nginx-logs
        pcaps_input=$input_dir/pcaps
        suffix=".min"
    fi
done

export BENCHMARK_CATEGORY="log-analysis"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

echo "nginx"
BENCHMARK_INPUT_FILE="$(realpath "$nginx_input")"
export BENCHMARK_INPUT_FILE
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/nginx.sh")"
export BENCHMARK_SCRIPT

$BENCHMARK_SHELL $scripts_dir/nginx.sh $nginx_input $outputs_dir/nginx$suffix 
echo $?
 
echo "pcaps"
BENCHMARK_INPUT_FILE="$(realpath "$pcaps_input")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/pcaps.sh")"
export BENCHMARK_SCRIPT

$BENCHMARK_SHELL $scripts_dir/pcaps.sh $pcaps_input $outputs_dir/pcaps$suffix 
echo $?
