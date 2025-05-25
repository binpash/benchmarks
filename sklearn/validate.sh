#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/sklearn"
# shell script to run verify.py
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
TMP="$eval_dir/inputs/input_$size"
export TMP
OUT="$eval_dir/outputs/out_$size"
export OUT
# run the Python script
python3 validate.py "${parsed_args[@]}"

# check if the script ran successfully
echo "sklearn $?"