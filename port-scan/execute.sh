#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/port-scan"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"

size="full"
for arg in "$@"; do
    case "$arg" in
        --small) size="small" ;;
        --min) size="min" ;;
    esac
done
export LC_ALL=C
inputs_dir="$inputs_dir/$size"
outputs_dir="$eval_dir/outputs/$size"
mkdir -p "$outputs_dir"

export BENCHMARK_CATEGORY="port-scan"
export BENCHMARK_SCRIPT="$scripts_dir/port-scan.sh"
export BENCHMARK_INPUT_FILE="$inputs_dir/port-scan.log"
go_install_dir="${eval_dir}/go_install"

export PATH=$PATH:/$go_install_dir/go/bin
export GOPATH=$HOME/go
export PATH=$PATH:$GOPATH/bin

$BENCHMARK_SHELL "$scripts_dir/port-scan.sh" "$inputs_dir/port-scan.log" "$eval_dir/inputs/routeviews.mrt" "$outputs_dir"

echo $?