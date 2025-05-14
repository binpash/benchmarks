#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done
REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/ray-tracing/outputs"
input_dir="${REPO_TOP}/ray-tracing/inputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi
