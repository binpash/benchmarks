#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${REPO_TOP}/git-workflow"
REPO_PATH="${EVAL_DIR}/inputs/chromium"
COMMITS_DIR="${EVAL_DIR}/inputs/commits"

cd "$REPO_PATH" || { echo "Cannot cd into $REPO_PATH"; exit 1; }

branch=$(git rev-parse --abbrev-ref HEAD)
if [ "$branch" != "bench_branch" ]; then
    echo "[ERROR] Expected to be on 'bench_branch', but found '$branch'"
    exit 1
fi

base_commit_file="$COMMITS_DIR/base_commit.txt"
if [ ! -f "$base_commit_file" ]; then
    echo "[ERROR] Missing base commit file at $base_commit_file"
    exit 1
fi

base_commit=$(cat "$base_commit_file")

commit_count=$(git rev-list --count "$base_commit"..HEAD)
if [ "$commit_count" -lt 5 ]; then
    echo "[ERROR] Expected at least 5 new commits after base commit, found $commit_count"
    exit 1
fi
