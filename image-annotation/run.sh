#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/image-annotation"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"
source venv/bin/activate

export BENCHMARK_CATEGORY="image-annotation"

export BENCHMARK_SCRIPT="$scripts_dir/image-annotation.sh"

suffix=""
for arg in "$@"; do
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done
inputs_dir="$inputs_dir/jpg$suffix"

BENCHMARK_INPUT_FILE="$(realpath "$inputs_dir")"
export BENCHMARK_INPUT_FILE

$BENCHMARK_SHELL "$scripts_dir/image-annotation.sh" "$inputs_dir" "$outputs_dir"

echo $?