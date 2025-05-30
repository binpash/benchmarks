#!/bin/bash
set -eu

BASE_DIR="$(dirname "$(readlink -f "$0")")"
status=0
if grep -q "FAIL" "${BASE_DIR}/run_results.log" 2>/dev/null; then
    status=1
    echo makeself $status
    exit $status
fi
echo makeself $status