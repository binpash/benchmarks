#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/sklearn"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="sklearn"
 
cd "$eval_dir" || exit 1 # scripts/execute.sh references PWD
 
BENCHMARK_SCRIPT="$(realpath "$scripts_dir/execute.sh")"
export BENCHMARK_SCRIPT
parsed_args=()
for arg in "$@"; do
    case "$arg" in
        --small|--min)
            parsed_args+=("$arg")
            ;;
    esac
done
$KOALA_SHELL "$scripts_dir/execute.sh" "${parsed_args[@]}"

