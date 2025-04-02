#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/opt-parallel"
scripts_dir="${eval_dir}/scripts"

SUITE_DIR="$(realpath "$(dirname "$0")")"
export SUITE_DIR
cd "$SUITE_DIR" || exit 1

export TIMEFORMAT=%R

suffix=""

for arg in "$@"; do
    case "$arg" in
        --small) suffix="_small" ;;
        --min) suffix="_min" ;;
    esac
done

export input_dir="${eval_dir}/inputs${suffix}/ChessData"

echo "executing opt-parallel $(date)"
mkdir -p "${SUITE_DIR}/outputs"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="opt-parallel"

main_script="${scripts_dir}/opt-parallel.sh"
BENCHMARK_SCRIPT="$(realpath "$main_script")"
export BENCHMARK_SCRIPT

"$BENCHMARK_SHELL" "$main_script" "$input_dir" > "${SUITE_DIR}/outputs/opt-parallel${suffix}.out"
