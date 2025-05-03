#!/bin/bash
set -eu

BASE_DIR="$(dirname "$(readlink -f "$0")")"

if grep -q "FAIL" "${BASE_DIR}/run_results.log" 2>/dev/null; then
    exit 1
fi