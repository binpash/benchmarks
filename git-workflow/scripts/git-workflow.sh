#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${REPO_TOP}/git-workflow"
REPO_PATH="${EVAL_DIR}/inputs/chromium"
COMMITS_DIR="${EVAL_DIR}/inputs/commits"
OUTPUT_DIR="${1:-${EVAL_DIR}/outputs}"

mkdir -p "$OUTPUT_DIR"
cd "$REPO_PATH" || { echo "Failed to enter repo at $REPO_PATH"; exit 1; }

git checkout benchmark-base
git reset --hard
git clean -fdx

for i in $(seq 5 -1 1); do
    statustime="/usr/bin/time --output=${OUTPUT_DIR}/status-${i}.log --verbose"
    addtime="/usr/bin/time --output=${OUTPUT_DIR}/add-${i}.log --verbose"
    committime="/usr/bin/time --output=${OUTPUT_DIR}/commit-${i}.log --verbose"
    patchtime="/usr/bin/time --output=${OUTPUT_DIR}/patchprep-${i}.log --verbose"

    patchfile="${COMMITS_DIR}/${i}-$((i-1)).diff"
    commitmsg="${COMMITS_DIR}/$((i-1)).commit"

    ${patchtime} patch -p1 < "$patchfile"

    ${statustime} git status

    ${addtime} git add .

    ${committime} git commit -F "$commitmsg"
done