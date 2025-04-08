#!/bin/bash
echo "Python: $(which python3)"
echo "Python version: $(python3 --version)"
echo "VIRTUAL_ENV: $VIRTUAL_ENV"

REPO_TOP="."
eval_dir="."
input_dir="${eval_dir}/input"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
mkdir -p "$results_dir"

img_input="${input_dir}/images_full"
suffix=".full"

if [[ "$@" == *"--small"* ]]; then
    img_input="${input_dir}/images_small"
    suffix=".small"
fi

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

echo "dpt sequential"
$BENCHMARK_SHELL "$scripts_dir/dpt_seq.sh" "$img_input" "$results_dir/seq_output$suffix.txt" \
    > "$results_dir/seq_log$suffix.txt"

echo "$?"

echo "dpt parallel"
$BENCHMARK_SHELL "$scripts_dir/dpt_par.sh" "$img_input" "$results_dir/par_output$suffix.txt" \
    > "$results_dir/par_log$suffix.txt"

echo "$?"