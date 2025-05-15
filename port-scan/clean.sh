#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/port-scan/outputs"
input_dir="${REPO_TOP}/port-scan/inputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
    rm -rf "${REPO_TOP}/port-scan/go_install"
fi
