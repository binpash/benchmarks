#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

REPO_TOP="$(git rev-parse --show-toplevel)"
input_dir="${REPO_TOP}/ci-cd/inputs/scripts/sqlite"
rm -rf "$input_dir"
