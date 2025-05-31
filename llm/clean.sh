#!/bin/bash
for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

TOP=$(git rev-parse --show-toplevel)
eval_dir="$TOP/llm"
input_dir="$TOP/llm/inputs"
outputs_dir="$eval_dir/outputs"

rm -rf "$outputs_dir"
rm -f "$eval_dir/ollama_serve.log"
if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi

