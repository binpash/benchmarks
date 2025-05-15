#!/bin/bash
for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done


REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/playlist-creation"
input_dir="$REPO_TOP/llm/inputs/scripts/playlist-creation/inputs"
outputs_dir="$eval_dir/outputs"
rm -rf "$outputs_dir"
if [ "$force" = true ]; then
    rm -rf "$input_dir"
    rm -rf "$eval_dir/venv"
fi
