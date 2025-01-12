#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/input"
scripts_dir="${eval_dir}/scripts"
results_dir="${eval_dir}/results"
mkdir -p $results_dir

nginx_input=$input_dir/nginx-logs
pcaps_input=$input_dir/pcaps
suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    # TODO: vary input size
    nginx_input=$input_dir/nginx-logs
    pcaps_input=$input_dir/pcaps
    suffix=".small"
fi

export BENCHMARK_CATEGORY="log-analysis"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
echo "shell: $BENCHMARK_SHELL"

echo "nginx"
export BENCHMARK_INPUT_FILE="$(realpath "$nginx_input")"
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/nginx.sh")"
time $BENCHMARK_SHELL $scripts_dir/nginx.sh $nginx_input $results_dir/nginx$suffix 
echo $?

echo "pcaps"
export BENCHMARK_INPUT_FILE="$(realpath "$pcaps_input")"
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/pcaps.sh")"
time $BENCHMARK_SHELL $scripts_dir/pcaps.sh $pcaps_input $results_dir/pcaps$suffix 
echo $?
