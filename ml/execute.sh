#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/ml"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="ml"
export LC_ALL=C

parsed_args=()
size="full"
for arg in "$@"; do
    case "$arg" in
        --small)
            parsed_args+=("$arg")
            size="small"
            ;;
        --min)
            parsed_args+=("$arg")
            size="min"
            ;;
    esac
done

echo "sklearn"
cd "$eval_dir" || exit 1 # scripts/execute.sh references PWD
 
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/execute.sh")"
export BENCHMARK_SCRIPT

TMP="$eval_dir/inputs/input_$size"
export TMP
OUT="$eval_dir/outputs/out_$size"
export OUT
mkdir -p "$OUT"
export BENCHMARK_INPUT_FILE="$eval_dir/inputs/input_$size"
$KOALA_SHELL "$scripts_dir/execute.sh" "${parsed_args[@]}"
echo "$?"

cd "$eval_dir" || exit 1

echo "dpt"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"

img_input="${input_dir}/images_full"
suffix=".full"

if [[ " $* " == *" --small "* ]]; then
    img_input="${input_dir}/images_small"
    suffix=".small"
fi
if [[ " $* " == *" --min "* ]]; then
    img_input="${input_dir}/images_min"
    suffix=".min"
fi
export BENCHMARK_INPUT_FILE="$img_input"
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/dpt_seq.sh")"
$KOALA_SHELL "$scripts_dir/dpt_seq.sh" "$img_input" "$outputs_dir/seq_output$suffix.txt"
echo "$?"
