#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${REPO_TOP}/git-workflow"
REPO_PATH="${EVAL_DIR}/inputs/chromium"
COMMITS_DIR="${EVAL_DIR}/inputs/commits"
OUTPUT_DIR="${1:-${EVAL_DIR}/outputs}"
NUM_COMMITS="${2:-20}"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$COMMITS_DIR"

cd "$REPO_PATH" || exit 1

git checkout main
git branch -D bench_branch 2>/dev/null || true
git checkout -b bench_branch
git reset --hard
git clean -fd

total_commits=$((NUM_COMMITS + 1))
commit_list=($(git rev-list --first-parent HEAD -n "$total_commits" | tac))

echo "${commit_list[0]}" > "$COMMITS_DIR/base_commit.txt"

for i in $(seq 1 "$NUM_COMMITS"); do
    prev="${commit_list[$((i - 1))]}"
    curr="${commit_list[$i]}"
    patchfile="$COMMITS_DIR/${i}-$((i - 1)).diff"
    commitmsg="$COMMITS_DIR/$((i - 1)).commit"

    git diff "$prev" "$curr" > "$patchfile"
    git log -1 --pretty=%B "$curr" > "$commitmsg"
done

git reset --hard "${commit_list[0]}"
git clean -fdx

git status

if [ -f "$COMMITS_DIR/base_commit.txt" ]; then
    git checkout bench_branch
    git reset --hard "$(cat "$COMMITS_DIR/base_commit.txt")"
else
    echo "Missing base_commit.txt"
    exit 1
fi

for i in $(seq "$NUM_COMMITS" -1 1); do
    patchfile="$COMMITS_DIR/${i}-$((i - 1)).diff"
    commitmsg="$COMMITS_DIR/$((i - 1)).commit"

    # Apply the patch
    patch -p1 < "$patchfile"

    git status

    git add -A
    git commit -F "$commitmsg"
done
