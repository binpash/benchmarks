#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
mkdir -p "$outputs_dir"

img_convert_input="$input_dir/jpg_full/jpg"
to_mp3_input="$input_dir/wav_full"
suffix=".full"

for arg in "$@"; do
    if [ "$arg" = "--min" ]; then
        img_convert_input="$input_dir/jpg_min/jpg"
        to_mp3_input="$input_dir/wav_min"
        suffix=".min"
        break
    elif [ "$arg" = "--small" ]; then
        img_convert_input="$input_dir/jpg_small/jpg"
        to_mp3_input="$input_dir/wav_small"
        suffix=".small"
        break
    fi
done

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="media-conv"
 
echo "img_convert"
BENCHMARK_INPUT_FILE="$(realpath "$img_convert_input")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/img_convert.sh")"
export BENCHMARK_SCRIPT

$KOALA_SHELL "$scripts_dir/img_convert.sh" "$img_convert_input" "$outputs_dir/img_convert$suffix" > "$outputs_dir/img_convert$suffix.log"
echo $?
 
echo "to_mp3"
BENCHMARK_INPUT_FILE="$(realpath "$to_mp3_input")"
export BENCHMARK_INPUT_FILE
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/to_mp3.sh")"
export BENCHMARK_SCRIPT

$KOALA_SHELL "$scripts_dir/to_mp3.sh" "$to_mp3_input" "$outputs_dir/to_mp3$suffix" > "$outputs_dir/to_mp3$suffix.log"
echo $?
