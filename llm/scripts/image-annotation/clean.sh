#!/bin/bash
KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/image-annotation"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
rm -r "$inputs_dir" "$outputs_dir"