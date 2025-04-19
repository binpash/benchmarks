#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/ray-tracing"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"

export BENCHMARK_CATEGORY="ray-tracing"

export BENCHMARK_SCRIPT="$scripts_dir/1.sh"
export BENCHMARK_INPUT_FILE="$inputs_dir/1.INFO"
$BENCHMARK_SHELL "$scripts_dir/1.sh" "$inputs_dir" "$outputs_dir"

export BENCHMARK_SCRIPT="$scripts_dir/2.sh"
$BENCHMARK_SHELL "$scripts_dir/2.sh" "$inputs_dir" "$outputs_dir"

export BENCHMARK_SCRIPT="$scripts_dir/3.sh"
$BENCHMARK_SHELL "$scripts_dir/3.sh" "$inputs_dir" "$outputs_dir"


export BENCHMARK_SCRIPT="$scripts_dir/4.sh"
$BENCHMARK_SHELL "$scripts_dir/4.sh" "$outputs_dir/rays.csv" > "$outputs_dir/query_pathid.log"

export BENCHMARK_SCRIPT="$scripts_dir/5.sh"
$BENCHMARK_SHELL "$scripts_dir/5.sh" "$outputs_dir/rays.csv" > "$outputs_dir/query_max_hop.log"

echo $?