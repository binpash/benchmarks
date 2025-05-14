#!/bin/bash

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/image-annotation"
scripts_dir="$eval_dir/scripts"
input_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"

export BENCHMARK_SCRIPT="$scripts_dir/image-annotation.sh"

suffix=""
for arg in "$@"; do
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done
input_dir="$input_dir/jpg$suffix"

BENCHMARK_INPUT_FILE="$(realpath "$input_dir")"
export BENCHMARK_INPUT_FILE

$KOALA_SHELL "$scripts_dir/image-annotation.sh" "$input_dir" "$outputs_dir"

echo $?