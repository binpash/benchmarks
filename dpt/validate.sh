#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"
mkdir -p "$hashes_dir"

suffix=".full"

generate=false
for arg in "$@"; do
    case $arg in
        --small)
            suffix=".small"
            ;;
        --min)
            suffix=".min"
            ;;
        --generate)
            generate=true
            ;;
    esac
done

# reference
if $generate; then
    python3 clean_output.py "$outputs_dir/seq_output$suffix.txt" "$outputs_dir/seq_output$suffix-cleaned.txt"
    seq_hash=$(shasum -a 256 "$outputs_dir/seq_output$suffix-cleaned.txt" | awk '{ print $1 }')
#    par_hash=$(shasum -a 256 "$outputs_dir/par_output$suffix-cleaned.txt" | awk '{ print $1 }')
    echo "$seq_hash" > "$hashes_dir/seq_output$suffix.txt"
#    echo "$par_hash" > "$hashes_dir/par_output$suffix.txt"
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
