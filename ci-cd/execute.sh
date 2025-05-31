#!/usr/bin/env bash
set -eu

BASE_DIR="$(dirname "$(readlink -f "$0")")"
TESTS_DIR="${BASE_DIR}/makeself/test"
LOGFILE="${BASE_DIR}/run_results.log"
KOALA_SHELL="${KOALA_SHELL:-bash}"
export BENCHMARK_CATEGORY="ci-cd"

echo "Starting test execution..." > "${LOGFILE}"

for test_script in "${TESTS_DIR}"/*/*.sh; do
    test_dir="$(dirname "${test_script}")"
    test_name="$(basename "${test_dir}")"
    test_log="${test_dir}/test_results.log"

    echo "Running test: ${test_name}" >> "${LOGFILE}"
    BENCHMARK_SCRIPT="$(realpath "$test_script")"
    export BENCHMARK_SCRIPT
    if $KOALA_SHELL "${test_script}" >> "${test_log}" 2>&1; then
        echo "PASS: ${test_name}" >> "${LOGFILE}"
    else
        echo "FAIL: ${test_name}" >> "${LOGFILE}"
    fi
done

echo "Test execution completed. Results in ${LOGFILE}"

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"

KOALA_SHELL=${KOALA_SHELL:-bash}

small_benchmark=(
    "xz-clang"
)


run_small=false

for arg in "$@"; do
    if [ "$arg" = "--min" ]; then
        run_small=true
        break
    fi
done

if [ "$run_small" = true ]; then
    for bench in "${small_benchmark[@]}"; do
        script_path="$eval_dir/$bench/execute.sh"
        if [ -x "$script_path" ]; then
            export BENCHMARK_SCRIPT="$script_path"
            $KOALA_SHELL $script_path "$@"
        else
            echo "Error: $script_path not found or not executable."
            exit 1
        fi
    done
    exit 0
fi

for bench in "$eval_dir"/*; do
    export BENCHMARK_SCRIPT="$bench/execute.sh"
    $KOALA_SHELL "$bench/execute.sh" "$@"
done