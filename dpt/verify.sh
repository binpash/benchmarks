#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"
mkdir -p "$hashes_dir"
source venv/bin/activate

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
    seq_hash=$(shasum -a 256 "$outputs_dir/seq_output$suffix-cleaned.txt" | awk '{ print $1 }')
#    par_hash=$(shasum -a 256 "$outputs_dir/par_output$suffix-cleaned.txt" | awk '{ print $1 }')
    echo "$seq_hash" > "$hashes_dir/seq_output$suffix.txt"
#    echo "$par_hash" > "$hashes_dir/par_output$suffix.txt"
    exit 0
fi

python clean_output.py "$outputs_dir/seq_output$suffix.txt" "$outputs_dir/seq_output$suffix-cleaned.txt"
seq_hash=$(shasum -a 256 "$outputs_dir/seq_output$suffix-cleaned.txt" | awk '{ print $1 }')
expected_sec_hash=$(cat "$hashes_dir/seq_output$suffix.txt")
status=0
if [[ "$seq_hash" != "$expected_sec_hash" ]]; then
    status=1
fi
#diff -q "$outputs_dir/seq_output$suffix-cleaned.txt" "$hashes_dir/seq_output$suffix.txt"
echo "seq $status"

# echo "Comparing par vs reference"
# python clean_output.py "$outputs_dir/par_output$suffix.txt" "$outputs_dir/par_output$suffix-cleaned.txt"

# par_hash=$(shasum -a 256 "$outputs_dir/par_output$suffix-cleaned.txt" | awk '{ print $1 }')
# expected_par_hash=$(cat "$hashes_dir/par_output$suffix.txt")
# status=0
# if [[ "$par_hash" != "$expected_par_hash" ]]; then
#     status=1
# fi
#diff -q "$outputs_dir/par_output$suffix-cleaned.txt" "$hashes_dir/par_output$suffix.txt"
# echo "par $status"

# This will result in differences due to processing time
# echo "Comparing seq vs par reference"
# diff -q "$outputs_dir/seq_output$suffix-cleaned.txt" "$outputs_dir/par_output$suffix-cleaned.txt"
# echo "seq vs par $?"