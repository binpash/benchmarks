#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
INPUT_DIR="${eval_dir}/inputs"
CHROMIUM_DIR="${INPUT_DIR}/chromium"

cd "$CHROMIUM_DIR" || exit 1
git stash --include-untracked >/dev/null 2>&1

git checkout benchmark-base
git branch -D benchmark-temp 2>/dev/null || true

git branch -D test_branch 2>/dev/null || true
cd "$eval_dir" || exit 1
rm -rf inputs outputs
