#!/usr/bin/env bash
set -eu

BASE_DIR="$(dirname "$(readlink -f "$0")")"
TESTS_DIR="${BASE_DIR}/makeself/test"
LOGFILE="${BASE_DIR}/run_results.log"
KOALA_SHELL="${KOALA_SHELL:-bash}"
export BENCHMARK_CATEGORY="makeself"

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
