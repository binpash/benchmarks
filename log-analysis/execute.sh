#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/inputs"
scripts_dir="${eval_dir}/scripts"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"

nginx_input=$input_dir/nginx-logs
pcaps_input=$input_dir/pcaps
size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done

nginx_input=$input_dir/nginx-logs_$size
pcaps_input=$input_dir/pcaps_$size

export BENCHMARK_CATEGORY="log-analysis"
KOALA_SHELL=${KOALA_SHELL:-bash}

echo "nginx"
BENCHMARK_INPUT_FILE="$(realpath "$nginx_input")"
export BENCHMARK_INPUT_FILE
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/nginx.sh")"
export BENCHMARK_SCRIPT

$KOALA_SHELL $scripts_dir/nginx.sh $nginx_input $outputs_dir/nginx_$size 
echo $?
 
echo "pcaps"
BENCHMARK_INPUT_FILE="$(realpath "$pcaps_input")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/pcaps.sh")"
export BENCHMARK_SCRIPT

$KOALA_SHELL $scripts_dir/pcaps.sh $pcaps_input $outputs_dir/pcaps_$size
echo $?
