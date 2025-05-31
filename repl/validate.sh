#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${REPO_TOP}/git-workflow"
REPO_PATH="${EVAL_DIR}/inputs/chromium"
COMMITS_DIR="${EVAL_DIR}/inputs/commits"

cd "$REPO_PATH" || { echo "Cannot cd into $REPO_PATH"; exit 1; }
NUM_COMMITS=21

for arg in "$@"; do
    case "$arg" in
        --min) NUM_COMMITS=2 ;;
        --small) NUM_COMMITS=6 ;;
    esac
done

check_commits=$((NUM_COMMITS-1))

branch=$(git rev-parse --abbrev-ref HEAD)
if [ "$branch" != "bench_branch" ]; then
    echo "Expected to be on 'bench_branch', but found '$branch'"
    exit 1
fi

base_commit_file="$COMMITS_DIR/base_commit.txt"
if [ ! -f "$base_commit_file" ]; then
    echo "Missing base commit file at $base_commit_file"
    exit 1
fi

base_commit=$(cat "$base_commit_file")

commit_count=$(git rev-list --count "$base_commit"..HEAD)
if [ "$commit_count" -lt $check_commits ]; then
    echo "Expected at least $check_commits new commits after base commit, found $commit_count"
    exit 1
fi
echo git-workflow 0

GENERATE=false

for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        GENERATE=true
        break
    fi
done

if [[ "$GENERATE" == true ]]; then
    python3 utils/validate.py --generate
    exit 0
else
    python3 utils/validate.py
    echo "vps-audit $?"
fi

