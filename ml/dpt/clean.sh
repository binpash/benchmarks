#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done
REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/dpt/outputs"
input_dir="${REPO_TOP}/dpt/inputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi

