#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/opt-parallel"
scripts_dir="${eval_dir}/scripts"

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

suffix=""

if [[ "$@" == *"--small"* ]]; then
    suffix="_small"
elif [[ "$@" == *"--min"* ]]; then
    suffix="_min"
fi

export input_dir="${eval_dir}/input${suffix}/ChessData"

echo "executing opt-parallel $(date)"
mkdir -p "${SUITE_DIR}/outputs"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="opt-parallel"

main_script="${scripts_dir}/opt-parallel.sh"
export BENCHMARK_SCRIPT="$(realpath "$main_script")"

"$BENCHMARK_SHELL" "$main_script" "$input_dir" > "${SUITE_DIR}/outputs/opt-parallel${suffix}.out"
