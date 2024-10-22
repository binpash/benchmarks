#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
results_dir="${eval_dir}/results"
correct_dir="${eval_dir}/correct-results"

if [[ "$@" == *"--generate"* ]]; then
    cp -r $results_dir/* "$correct_dir"
    exit 0
fi

diff -q "$results_dir/average.txt" "$correct_dir/average.txt"
echo average $?

diff -q "$results_dir/min.txt" "$correct_dir/min.txt"
echo min $?

diff -q "$results_dir/max.txt" "$correct_dir/max.txt"
echo max $?

