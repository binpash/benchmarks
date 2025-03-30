#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
results_dir="${eval_dir}/results"
hashes_dir="${eval_dir}/hashes"
mkdir -p "$hashes_dir"

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

# reference
if [[ "$@" == *"--generate"* ]]; then
    cp "$results_dir/seq_output$suffix.txt" "$hashes_dir/seq_output$suffix.txt"
    cp "$results_dir/par_output$suffix.txt" "$hashes_dir/par_output$suffix.txt"
    echo "Done."
    exit 0
fi

echo "Comparing seq vs reference"
diff -q "$results_dir/seq_output$suffix.txt" "$hashes_dir/seq_output$suffix.txt"
echo "seq $?"

echo "Comparing par vs reference"
diff -q "$results_dir/par_output$suffix.txt" "$hashes_dir/par_output$suffix.txt"
echo "par $?"

echo "Comparing seq vs par reference"
diff -q "$results_dir/seq_output$suffix.txt" "$results_dir/par_output$suffix.txt"
echo "seq vs par $?"
