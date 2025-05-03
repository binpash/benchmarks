#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"

small_benchmark=(
    "lua"
    "memcached"
    "redis"
    "sqlite"
    "vim"
    "xz"
    "xz-clang"
)

run_small=false

for arg in "$@"; do
    if [ "$arg" = "--small" ]; then
        run_small=true
        break
    fi
done


if [ "$run_small" = true ]; then
    for bench in "${small_benchmark[@]}"; do
        script_path="$scripts_dir/$bench/fetch.sh"
        if [ -x "$script_path" ]; then
            "$script_path" "$@"
        else
            echo "Error: $script_path not found or not executable."
            exit 1
        fi
    done
    exit 0
fi

for bench in "$scripts_dir"/*; do
    "$bench/fetch.sh" "$@"
done

