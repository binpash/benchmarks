#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs/scripts/redis"
rm -rf "$input_dir/dev"
