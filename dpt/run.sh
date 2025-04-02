#!/bin/bash
source .venv/bin/activate

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo "VIRTUAL_ENV: $VIRTUAL_ENV"

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
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
