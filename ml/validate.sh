#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/ml"
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"
# shell script to run verify.py
parsed_args=()

suffix=".full"

size="full"
generate=false
for arg in "$@"; do
    case "$arg" in
        --small)
            parsed_args+=("$arg")
            size="small"
            suffix=".small"
            ;;
        --min)
            parsed_args+=("$arg")
            size="min"
            suffix=".min"
            ;;
        --generate)
            generate=true
            ;;
    esac
done
TMP="$eval_dir/inputs/input_$size"
export TMP
OUT="$eval_dir/outputs/out_$size"
export OUT

# run the Python script
python3 validate.py "${parsed_args[@]}"
echo "sklearn $?"

if $generate; then
    python3 clean_output.py "$outputs_dir/seq_output$suffix.txt" "$outputs_dir/seq_output$suffix-cleaned.txt"
    seq_hash=$(shasum -a 256 "$outputs_dir/seq_output$suffix-cleaned.txt" | awk '{ print $1 }')
    echo "$seq_hash" > "$hashes_dir/seq_output$suffix.txt"
    exit 0
fi
python3 clean_output.py "$outputs_dir/seq_output$suffix.txt" "$outputs_dir/seq_output$suffix-cleaned.txt"
seq_hash=$(shasum -a 256 "$outputs_dir/seq_output$suffix-cleaned.txt" | awk '{ print $1 }')
expected_sec_hash=$(cat "$hashes_dir/seq_output$suffix.txt")

status=0
if [[ "$seq_hash" != "$expected_sec_hash" ]]; then
    status=1
fi
echo "seq $status"