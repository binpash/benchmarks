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

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

echo "nginx"
$BENCHMARK_SHELL $scripts_dir/nginx.sh $nginx_input $results_dir/nginx$suffix 
echo $?

echo "pcaps"
$BENCHMARK_SHELL $scripts_dir/pcaps.sh $pcaps_input $results_dir/pcaps$suffix 
echo $?
