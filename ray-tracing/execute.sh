#!/bin/bash

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/ray-tracing"
scripts_dir="$eval_dir/scripts"
input_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"

size="full"
for arg in "$@"; do
    case "$arg" in
        --small) size="small" ;;
        --min) size="min" ;;
    esac
done
input_dir="$eval_dir/inputs/$size"
outputs_dir="$eval_dir/outputs/$size"
mkdir -p "$outputs_dir"

export BENCHMARK_CATEGORY="ray-tracing"

export BENCHMARK_SCRIPT="$scripts_dir/1.sh"
export BENCHMARK_INPUT_FILE="$input_dir/1.INFO"
$KOALA_SHELL "$scripts_dir/1.sh" "$input_dir" "$outputs_dir"

export BENCHMARK_INPUT_FILE="$input_dir"

export BENCHMARK_SCRIPT="$scripts_dir/2.sh"
$KOALA_SHELL "$scripts_dir/2.sh" "$input_dir" "$outputs_dir"

export BENCHMARK_SCRIPT="$scripts_dir/3.sh"
$KOALA_SHELL "$scripts_dir/3.sh" "$input_dir" "$outputs_dir"


export BENCHMARK_SCRIPT="$scripts_dir/4.sh"
$KOALA_SHELL "$scripts_dir/4.sh" "$outputs_dir/rays.csv" > "$outputs_dir/query_pathid.log"

export BENCHMARK_SCRIPT="$scripts_dir/5.sh"
$KOALA_SHELL "$scripts_dir/5.sh" "$outputs_dir/rays.csv" > "$outputs_dir/query_max_hop.log"

echo $?