#!/bin/bash

while getopts "f" opt; do
  case $opt in
    f) force=true ;;
  esac
done

REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/media-conv/outputs"
input_dir="${REPO_TOP}/media-conv/inputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi
