#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN="$REPO_TOP/aurpkg/inputs/packages"
OUT="${OUT:-$REPO_TOP/aurpkg/outputs}"
BENCHMARK_SHELL="${BENCHMARK_SHELL:-bash}"
SCRIPT="./scripts/pacaur.sh"

test "$UID" -gt 0 || { echo "Don't run this as root!"; exit 1; } 

mkdir -p "${OUT}"

# Set environment variables
BENCHMARK_CATEGORY="aurpkg"
export BENCHMARK_CATEGORY

BENCHMARK_SCRIPT="$(realpath "$SCRIPT")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN")"
export BENCHMARK_INPUT_FILE

echo "$SCRIPT"
$BENCHMARK_SHELL "$SCRIPT" "$IN" "$OUT"
echo "$?"
