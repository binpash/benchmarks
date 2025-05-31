#! /bin/bash

#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"

small_benchmark=(
    "xz-clang"
)

run_small=false

for arg in "$@"; do
    if [ "$arg" = "--small" ]; then
        run_small=true
        break
    fi
    if [ "$arg" = "--min" ]; then
        run_small=true
        break
    fi
done

if [ "$run_small" = true ]; then
    for bench in "${small_benchmark[@]}"; do
        script_path="$eval_dir/$bench/fetch.sh"
        if [ -x "$script_path" ]; then
            "$script_path" "$@"
        else
            echo "Error: $script_path not found or not executable."
            exit 1
        fi
    done
    exit 0
fi

for bench in "$eval_dir"/*; do
    "$bench/fetch.sh" "$@"
done