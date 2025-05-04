#!/bin/bash
KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
INPUT_DIR="${eval_dir}/inputs"
OUTPUT_DIR="${eval_dir}/outputs"
echo "INPUT_DIR: $INPUT_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
CHROMIUM_DIR="${INPUT_DIR}/chromium"

cd "$CHROMIUM_DIR" || exit 1
git stash
git checkout main
git branch -D bench_branch 2>/dev/null || true

rm -rf "$INPUT_DIR"
rm -rf "$OUTPUT_DIR"
