#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/port-scan"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"

export BENCHMARK_CATEGORY="port-scan"
export BENCHMARK_SCRIPT="$scripts_dir/port-scan.sh"
export BENCHMARK_INPUT_FILE="$inputs_dir/port-scan.log"
go_install_dir="${eval_dir}/go_install"

export PATH=$PATH:/$go_install_dir/go/bin
export GOPATH=$HOME/go
export PATH=$PATH:$GOPATH/bin

$BENCHMARK_SHELL "$scripts_dir/port-scan.sh" "$inputs_dir/port-scan.log" "$inputs_dir/routeviews.mrt" "$outputs_dir"

echo $?
