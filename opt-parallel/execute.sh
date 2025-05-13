#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/opt-parallel"
scripts_dir="${eval_dir}/scripts"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"

size=full
input_dir="${eval_dir}/inputs/ChessData"
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done

export LC_ALL=C
if [ "$size" = "min" ]; then
    input_dir="${eval_dir}/inputs/ChessData_min"
fi

if [ "$size" = "small" ]; then
    input_dir="${eval_dir}/inputs/ChessData_small"
fi

export BENCHMARK_INPUT_FILE="${input_dir}"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="opt-parallel"

main_script="${scripts_dir}/opt-parallel.sh"
BENCHMARK_SCRIPT="$(realpath "$main_script")"
export BENCHMARK_SCRIPT

$KOALA_SHELL "$main_script" "$input_dir" > "${outputs_dir}/opt-parallel_$size.out"
