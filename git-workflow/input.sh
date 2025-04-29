#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
INPUT_DIR="${eval_dir}/inputs"
CHROMIUM_DIR="${INPUT_DIR}/chromium"
COMMITS_DIR="${INPUT_DIR}/commits"

mkdir -p "$CHROMIUM_DIR" "$COMMITS_DIR"

NUM_COMMITS=21

for arg in "$@"; do
    case "$arg" in
        --min) NUM_COMMITS=2 ;;
        --small) NUM_COMMITS=6 ;;
    esac
done

if [ ! -d "$CHROMIUM_DIR/.git" ]; then
    git clone --depth=$((NUM_COMMITS+2)) https://chromium.googlesource.com/chromium/src.git "$CHROMIUM_DIR" || {
        echo "Failed to clone Chromium repository"
        exit 1
    }
fi
