#! /bin/bash

BASE_DIR="$(dirname "$(readlink -f "$0")")"
TEST_DIR="${BASE_DIR}/makeself/test"
find "${TEST_DIR}" -type f -name "*.log" -exec rm -f {} +

if [[ -f run_results.log ]]; then
    rm run_results.log
fi

if [[ -f verify_results.log ]]; then
    rm verify_results.log
fi