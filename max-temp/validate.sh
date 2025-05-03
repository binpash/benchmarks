#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"

suffix=".full"

generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
    --small) suffix=".small" ;;
    --min) suffix=".min" ;;
    esac
done

statistics_dir="${eval_dir}/outputs/statistics$suffix"
correct_dir="${eval_dir}/correct-results/statistics$suffix"

if $generate; then
    mkdir -p "$correct_dir"
    cp -r "$statistics_dir"/* "$correct_dir"
    exit 0
fi

diff -q "$statistics_dir/average.txt" "$correct_dir/average.txt"
echo average$suffix $?

diff -q "$statistics_dir/min.txt" "$correct_dir/min.txt"
echo min$suffix $?

diff -q "$statistics_dir/max.txt" "$correct_dir/max.txt"
echo max$suffix $?
