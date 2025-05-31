#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

REPO_TOP=$(git rev-parse --show-toplevel)
input_dir="${REPO_TOP}/analytics/inputs"
outputs_dir="${REPO_TOP}/analytics/outputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
    rm -rf "${REPO_TOP}/analytics/go_install"
fi
