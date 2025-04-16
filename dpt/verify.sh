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
    cp "$outputs_dir/seq_output$suffix.txt" "$hashes_dir/seq_output$suffix.txt"
    cp "$outputs_dir/par_output$suffix.txt" "$hashes_dir/par_output$suffix.txt"
    echo "Done."
    exit 0
fi

echo "Comparing seq vs reference"
diff -q "$outputs_dir/seq_output$suffix.txt" "$hashes_dir/seq_output$suffix.txt"
echo "seq $?"

echo "Comparing par vs reference"
diff -q "$outputs_dir/par_output$suffix.txt" "$hashes_dir/par_output$suffix.txt"
echo "par $?"

# This will result in differences due to processing time
# echo "Comparing seq vs par reference"
# diff -q "$outputs_dir/seq_output$suffix.txt" "$outputs_dir/par_output$suffix.txt"
# echo "seq vs par $?"