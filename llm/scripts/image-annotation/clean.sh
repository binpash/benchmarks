#!/bin/bash
for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/image-annotation"
input_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi
