#!/bin/bash
set -eu

BASE_DIR="$(dirname "$(readlink -f "$0")")"
TESTS_DIR="${BASE_DIR}/makeself/test"

all_passed=true

for test_log in "${TESTS_DIR}"/*/test_results.log; do
    if grep -q "FAIL" "${test_log}" 2>/dev/null; then
        all_passed=false
        break
    fi
done

if [ "$all_passed" = true ]; then
    exit 0
else
    exit 1
fi
