#!/bin/bash
# source: posh benchmark suite

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

commit_file="$COMMITS_DIR/commit_list.txt"
git rev-list --first-parent HEAD -n "$NUM_COMMITS" | tac > "$commit_file"

base_commit=$(head -n 1 "$commit_file")
echo "$base_commit" > "$COMMITS_DIR/base_commit.txt"

num_patches=$((NUM_COMMITS - 1))
read -r base_commit < "$commit_file"
prev_commit="$base_commit"
i=1

while read -r curr_commit; do
    patch_upper=$(( num_patches - i + 1 ))
    patch_lower=$(( patch_upper - 1 ))
    
    patchfile="$COMMITS_DIR/${patch_upper}-${patch_lower}.diff"
    commitmsg="$COMMITS_DIR/${patch_upper}-${patch_lower}.commit"
    
    git diff "$prev_commit" "$curr_commit" > "$patchfile"
    git log -1 --pretty=%B "$curr_commit" > "$commitmsg"
    
    prev_commit="$curr_commit"
    i=$((i + 1))
done < <(tail -n +2 "$commit_file")


git add -A

git status

if [ -f "$COMMITS_DIR/base_commit.txt" ]; then
    git checkout bench_branch
    git reset --hard "$(cat "$COMMITS_DIR/base_commit.txt")"
else
    echo "Missing base_commit.txt"
    exit 1
fi

for i in $(seq "$num_patches" -1 1); do
    lower=$(( i - 1 ))
    patchfile="$COMMITS_DIR/${i}-${lower}.diff"
    commitmsg="$COMMITS_DIR/${i}-${lower}.commit"
    
    patch -p1 < "$patchfile" || { echo "Failed to apply $patchfile"; exit 1; }
    git status

    git add -A
    git commit -F "$commitmsg" || { echo "Failed to commit with $commitmsg"; exit 1; }
done


# Final status check
git status