#!/bin/bash
set -eu

BASE_DIR="$(dirname "$(readlink -f "$0")")"
TESTS_DIR="${BASE_DIR}/makeself/test"
VERIFY_LOG="${BASE_DIR}/verify_results.log"

echo "Verifying test results..." > "${VERIFY_LOG}"

all_passed=true

for test_log in "${TESTS_DIR}"/*/test_results.log; do
    test_name="$(basename "$(dirname "${test_log}")")"

    echo "Verifying: ${test_name}" >> "${VERIFY_LOG}"
    if grep -q "FAIL" "${test_log}"; then
        all_passed=false
        echo "FAIL: ${test_name}" >> "${VERIFY_LOG}"
        grep "FAIL" "${test_log}" >> "${VERIFY_LOG}"
    else
        echo "PASS: ${test_name}" >> "${VERIFY_LOG}"
    fi
done

if [[ "${all_passed}" == "true" ]]; then
    echo "All tests passed successfully!" >> "${VERIFY_LOG}"
    echo "All tests passed successfully!"
    exit 0
else
    echo "Some tests failed. Review the logs above." >> "${VERIFY_LOG}"
    echo "Some tests failed. Review the logs above."
    exit 1
fi
