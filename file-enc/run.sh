#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
mkdir -p "$outputs_dir"

echo "Executing file-enc benchmark at $(date)"

input_pcaps="$input_dir/pcaps"
suffix=".full"
for arg in "$@"; do
    case "$arg" in
        --small)
            suffix=".small"
            ;;
        --min)
            suffix=".min"
            ;;
    esac
done

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

BENCHMARK_CATEGORY="file-enc"
export BENCHMARK_CATEGORY

BENCHMARK_INPUT_FILE="$(realpath "$input_pcaps")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/compress_files.sh")"
export BENCHMARK_SCRIPT
$BENCHMARK_SHELL "$scripts_dir/compress_files.sh" "$input_pcaps" "$outputs_dir/compress_files$suffix"

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/encrypt_files.sh")"
export BENCHMARK_SCRIPT
$BENCHMARK_SHELL "$scripts_dir/encrypt_files.sh" "$input_pcaps" "$outputs_dir/encrypt_files$suffix"