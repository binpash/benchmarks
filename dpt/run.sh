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
export LC_ALL=C
export BENCHMARK_CATEGORY="dpt"
img_input="${input_dir}/images_full"
suffix=".full"

if [[ " $* " == *" --small "* ]]; then
    img_input="${input_dir}/images_small"
    suffix=".small"
fi

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_INPUT_FILE="$img_input"
echo "dpt sequential"
$BENCHMARK_SHELL "$scripts_dir/dpt_seq.sh" "$img_input" "$outputs_dir/seq_output$suffix.txt"
echo "$?"
cat "$outputs_dir/seq_output$suffix.txt"
echo "$(cat "$outputs_dir/seq_output$suffix.txt")"
# echo "dpt parallel"
# $BENCHMARK_SHELL "$scripts_dir/dpt_par.sh" "$img_input" "$outputs_dir/par_output$suffix.txt"
# echo "$?"