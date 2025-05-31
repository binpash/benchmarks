#!/bin/bash

cd "$(dirname "$0")" || exit 1

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

eval_dir="$PWD"
outputs_dir="${eval_dir}/outputs"
input_dir="${eval_dir}/inputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi

