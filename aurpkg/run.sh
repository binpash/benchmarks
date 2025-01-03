#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
IN=$REPO_TOP/aurpkg/input/packages
OUT=${OUT:-$REPO_TOP/aurpkg/outputs}
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
mkdir -p ${OUT}

script="./scripts/pacaur.sh"

# Switch to user "user" to avoid permission issues

echo "$script"
sudo -u user $BENCHMARK_SHELL "$script" "$IN" "$OUT"
echo "$?"
