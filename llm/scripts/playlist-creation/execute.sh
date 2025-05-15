#!/bin/bash

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/playlist-creation"
scripts_dir="$eval_dir/scripts"
input_dir="$REPO_TOP/llm/inputs/scripts/playlist-creation/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"

export BENCHMARK_SCRIPT="$scripts_dir/playlist-creation.sh"
source "$eval_dir/venv/bin/activate"

suffix=""
for arg in "$@"; do
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done
export LC_ALL=C
input_dir="$input_dir/songs$suffix"
outputs_dir="$outputs_dir/songs$suffix"
mkdir -p "$outputs_dir"
echo "input_dir: $input_dir"
BENCHMARK_INPUT_FILE="$(realpath "$input_dir")"
export BENCHMARK_INPUT_FILE

$KOALA_SHELL "$scripts_dir/playlist-creation.sh" "$input_dir" "$outputs_dir"

echo $?