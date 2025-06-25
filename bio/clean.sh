#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/bio"
outputs_dir="${eval_dir}/outputs"
input_dir="${eval_dir}/inputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi
