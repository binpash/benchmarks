#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/image_annotation"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"

export BENCHMARK_CATEGORY="image-annotation"

export BENCHMARK_SCRIPT="$scripts_dir/image_annotation.sh"

$BENCHMARK_SHELL "$scripts_dir/image_annotation.sh" "$inputs_dir" "$outputs_dir"

echo $?