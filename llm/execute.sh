#!/bin/bash

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm"
scripts_dir="$eval_dir/scripts"
input_dir="$REPO_TOP/llm/inputs"
outputs_dir="$eval_dir/outputs"

suffix=""
for arg in "$@"; do
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done
export LC_ALL=C

export BENCHMARK_CATEGORY="llm"
img_input_dir="$input_dir/jpg$suffix"
img_outputs_dir="$outputs_dir/jpg$suffix"
mkdir -p "$img_outputs_dir"

BENCHMARK_INPUT_FILE="$(realpath "$img_input_dir")"
export BENCHMARK_INPUT_FILE

export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/image-annotation.sh")"
$KOALA_SHELL "$scripts_dir/image-annotation.sh" "$img_input_dir" "$img_outputs_dir"
echo $?

songs_input_dir="$input_dir/songs$suffix"
songs_outputs_dir="$outputs_dir/songs$suffix"
mkdir -p "$songs_outputs_dir"

BENCHMARK_INPUT_FILE="$(realpath "$songs_input_dir")"
export BENCHMARK_INPUT_FILE

export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/playlist-creation.sh")"
$KOALA_SHELL "$scripts_dir/playlist-creation.sh" "$songs_input_dir" "$songs_outputs_dir"
echo $?