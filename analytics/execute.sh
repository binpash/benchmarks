#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/inputs"
scripts_dir="${eval_dir}/scripts"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"

export LC_ALL=C

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
port_scan_input=$input_dir/port_scan_$size/all_logs.jsonl
rt_inputs_dir=$input_dir/ray_tracing_$size
rt_outputs_dir=$outputs_dir/ray_tracing_$size

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

echo "port-scan"
BENCHMARK_INPUT_FILE="$(realpath "$port_scan_input")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/port-scan.sh")"
export BENCHMARK_SCRIPT
go_install_dir="${eval_dir}/go_install"

export PATH=$PATH:/$go_install_dir/go/bin
export GOPATH=$HOME/go
export PATH=$PATH:$GOPATH/bin

mkdir -p "$outputs_dir/port_scan_$size"
touch "$outputs_dir/port_scan_$size/annotated.jsonl"
touch "$outputs_dir/port_scan_$size/file1"
touch "$outputs_dir/port_scan_$size/file2"
touch "$outputs_dir/port_scan_$size/as_popularity.csv"

$KOALA_SHELL $scripts_dir/port-scan.sh "$port_scan_input" "$eval_dir/inputs/routeviews.mrt" "$outputs_dir/port_scan_$size/annotated" "$outputs_dir/port_scan_$size/file1" "$outputs_dir/port_scan_$size/file2" "$outputs_dir/port_scan_$size/as_popularity"
echo $?

echo "ray-tracing"
mkdir -p "$outputs_dir/ray_tracing_$size"
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/ray-tracing.sh")"
export BENCHMARK_SCRIPT
BENCHMARK_INPUT_FILE="$(realpath "$rt_inputs_dir")"
export BENCHMARK_INPUT_FILE
$KOALA_SHELL "$scripts_dir/ray-tracing.sh" "$rt_inputs_dir" "$rt_outputs_dir"
echo $?