#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
INPUT_DIR="${eval_dir}/inputs"
CHROMIUM_DIR="${INPUT_DIR}/chromium"
PATCH_DIR="${INPUT_DIR}/commits"

mkdir -p "$CHROMIUM_DIR" "$PATCH_DIR"

# Clone Chromium repo (shallow) if not already present
if [ ! -d "$CHROMIUM_DIR/.git" ]; then
    git clone --depth=30 https://chromium.googlesource.com/chromium/src.git "$CHROMIUM_DIR"
fi

cd "$CHROMIUM_DIR" || exit 1

if git show-ref --verify --quiet refs/heads/benchmark-base; then
    git checkout benchmark-base
else
    git checkout -b benchmark-base
fi

git checkout -B benchmark-temp benchmark-base

BENCH_FILE="benchmark.txt"
if [ ! -f "$BENCH_FILE" ]; then
    echo "Benchmark file initialized." > "$BENCH_FILE"
    git add "$BENCH_FILE"
    git commit -m "Initialize benchmark file"
fi

prev_commit=$(git rev-parse HEAD)


for i in {1..5}; do
    echo "Benchmark commit ${i}: $(date)" >> "$BENCH_FILE"
    
    git add "$BENCH_FILE"
    commit_msg="Benchmark commit ${i}"
    git commit -m "$commit_msg"
    
    cur_commit=$(git rev-parse HEAD)
    
    patch_file="${PATCH_DIR}/${i}-${i-1}.diff"
    git diff "$prev_commit" "$cur_commit" > "$patch_file"
    
    commit_file="${PATCH_DIR}/$((i-1)).commit"
    git log -1 --pretty=%B "$cur_commit" > "$commit_file"
    
    prev_commit="$cur_commit"
done

git checkout benchmark-base
