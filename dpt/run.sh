#!/bin/bash
echo "Python: $(which python3)"
echo "Python version: $(python3 --version)"
echo "VIRTUAL_ENV: $VIRTUAL_ENV"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/dpt"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
mkdir -p "$outputs_dir"
source venv/bin/activate

export BENCHMARK_CATEGORY="dpt"
export BENCHMARK_INPUT_FILE="$input_dir"
img_input="${input_dir}/images_full"
suffix=".full"

if [[ " $* " == *" --small "* ]]; then
    img_input="${input_dir}/images_small"
    suffix=".small"
fi

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

echo "dpt sequential"
$BENCHMARK_SHELL "$scripts_dir/dpt_seq.sh" "$img_input" "$outputs_dir/seq_output$suffix.txt" \
    > "$outputs_dir/seq_log$suffix.txt"

echo "$?"

echo "dpt parallel"
$BENCHMARK_SHELL "$scripts_dir/dpt_par.sh" "$img_input" "$outputs_dir/par_output$suffix.txt" \
    > "$outputs_dir/par_log$suffix.txt"

echo "$?"