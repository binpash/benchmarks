#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
INPUT_DIR="${eval_dir}/inputs"
CHROMIUM_DIR="${INPUT_DIR}/chromium"
COMMITS_DIR="${INPUT_DIR}/commits"

mkdir -p "$CHROMIUM_DIR" "$COMMITS_DIR"

if [ ! -d "$CHROMIUM_DIR/.git" ]; then
    git clone https://chromium.googlesource.com/chromium/src.git "$CHROMIUM_DIR" || {
        echo "Failed to clone Chromium repository"
        exit 1
    }
fi
