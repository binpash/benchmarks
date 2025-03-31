#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN="$REPO_TOP/aurpkg/inputs/packages"
OUT="${OUT:-$REPO_TOP/aurpkg/outputs}"
BENCHMARK_SHELL="${BENCHMARK_SHELL:-bash}"
SCRIPT="./scripts/pacaur.sh"

mkdir -p "$OUT"

# Switch to user "user" to avoid permission issues

# Set environment variables
BENCHMARK_CATEGORY="aurpkg"
export BENCHMARK_CATEGORY

BENCHMARK_SCRIPT="$(realpath "$SCRIPT")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN")"
export BENCHMARK_INPUT_FILE

echo "$SCRIPT"
"$BENCHMARK_SHELL" "$SCRIPT" "$IN" "$OUT"
echo "$?"
