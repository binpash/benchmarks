#!/bin/bash
for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/image-annotation"
input_dir="$REPO_TOP/llm/inputs/scripts/image-annotation/inputs"
outputs_dir="$eval_dir/outputs"

rm -rf "$outputs_dir"
rm -f "$eval_dir/ollama_serve.log"
if [ "$force" = true ]; then
    rm -rf "$input_dir"
    rm -rf "$eval_dir/venv"
fi
