#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd"

min_benchmark=(
    "xz-clang"
)

run_min=false

for arg in "$@"; do
    if [ "$arg" = "--min" ]; then
        run_min=true
        break
    fi
done

if [ "$run_min" = true ]; then
    for bench in "${min_benchmark[@]}"; do
        script_path="$eval_dir/riker/$bench/validate.sh"
        if [ -x "$script_path" ]; then
            "$script_path" "$@"
        else
            echo "Error: $script_path not found or not executable."
            exit 1
        fi
    done
else
    for bench in "$eval_dir"/riker/*; do
        "$bench/validate.sh" "$@"
    done
fi

status=0
if grep -q "FAIL" "${eval_dir}/run_results.log" 2>/dev/null; then
    status=1
    echo makeself $status
    exit $status
fi
echo makeself $status