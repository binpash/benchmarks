#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
results_dir="${eval_dir}/results"
correct_dir="${eval_dir}/correct-results"

diff "$results_dir/average.txt" "$correct_dir/average.txt" \
    && diff "$results_dir/min.txt" "$correct_dir/min.txt" \
    && diff "$results_dir/max.txt" "$correct_dir/max.txt"

if [ $? -eq 0 ]; then
    echo "Valid"
else
    echo "Invalid"
    exit 1
fi
