#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/file-mod"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
mkdir -p "$outputs_dir"

img_convert_input="$input_dir/jpg_full/jpg"
to_mp3_input="$input_dir/wav_full"
suffix=".full"

size="full"
input_pcaps="$input_dir/pcaps"
for arg in "$@"; do
    case "$arg" in
        --small)
            img_convert_input="$input_dir/jpg_small/jpg"
            to_mp3_input="$input_dir/wav_small"
            suffix=".small"
            size="small"
            ;;
        --min)
            size="small"
            img_convert_input="$input_dir/jpg_min/jpg"
            to_mp3_input="$input_dir/wav_min"
            suffix=".min"
            size="min"
            ;;
    esac
done

input_pcaps="$input_dir/pcaps_$size"


KOALA_SHELL=${KOALA_SHELL:-bash}

BENCHMARK_CATEGORY="file-mod"
export BENCHMARK_CATEGORY

BENCHMARK_INPUT_FILE="$(realpath "$input_pcaps")"
export BENCHMARK_INPUT_FILE

echo "compress_files"
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/compress_files.sh")"
export BENCHMARK_SCRIPT
$KOALA_SHELL "$scripts_dir/compress_files.sh" "$input_pcaps" "$outputs_dir/compress_files$suffix"
echo $?

echo "encrypt_files"
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/encrypt_files.sh")"
export BENCHMARK_SCRIPT
$KOALA_SHELL "$scripts_dir/encrypt_files.sh" "$input_pcaps" "$outputs_dir/encrypt_files$suffix"
echo $?

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

echo "thumbnail_generation"
BENCHMARK_INPUT_FILE="$(realpath "$img_convert_input")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "$scripts_dir/thumbnail_generation.sh")"
export BENCHMARK_SCRIPT
mkdir -p "$outputs_dir/thumbnail$suffix"
$KOALA_SHELL "$scripts_dir/thumbnail_generation.sh" "$img_convert_input" "$outputs_dir/thumbnail$suffix" > "$outputs_dir/thumbnail$suffix.log"
echo $?
