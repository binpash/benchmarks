#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/sklearn"
size="full"
parsed_args=()
size="full"
for arg in "$@"; do
    case "$arg" in
        --small)
            parsed_args+=("$arg")
            size="small"
            ;;
        --min)
            parsed_args+=("$arg")
            size="min"
            ;;
    esac
done

mkdir -p $eval_dir/outputs
mkdir -p $eval_dir/inputs/input_"$size"

export TMP="$eval_dir/inputs/input_$size"
# Generating model & samples
python3 $eval_dir/scripts/gen_model.py 100
python3 $eval_dir/scripts/gen_samples.py "${parsed_args[@]}"