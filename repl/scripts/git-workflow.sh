#!/bin/bash
# source: posh benchmark suite

shopt -s expand_aliases

alias g='git'
alias gst='git status'
alias gco='git checkout'
alias grs='git reset --hard'
alias gcl='git clean -fd'
alias gci='git commit'
alias gaa='git add -A'

TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${TOP}/repl"
REPO_PATH="${EVAL_DIR}/inputs/chromium"
COMMITS_DIR="${EVAL_DIR}/inputs/commits"
NUM_COMMITS=${NUM_COMMITS:-2}

# override HOME variable
export HOME="$COMMITS_DIR"
mkdir -p "$COMMITS_DIR"

cd "$REPO_PATH" || exit 1
git config --global --add safe.directory "$REPO_PATH"
g config user.email "author@example.com"
g config user.name "A U Thor"

g stash
gco main
g branch -D bench_branch 2>/dev/null || true
gco -b bench_branch
grs
gcl

commit_file=~/commit_list.txt
g rev-list --first-parent HEAD -n "$NUM_COMMITS" | tac > "$commit_file"

base_commit=$(head -n 1 "$commit_file")
echo "$base_commit" > ~/base_commit.txt

num_patches=$((NUM_COMMITS - 1))
read -r base_commit < "$commit_file"
prev_commit="$base_commit"
i=1

tail -n +2 "$commit_file" | while read -r curr_commit; do
    patch_upper=$((num_patches - i + 1))
    patch_lower=$((patch_upper - 1))
    
    patchfile=~/${patch_upper}-${patch_lower}.diff
    commitmsg=~/${patch_upper}-${patch_lower}.commit
    
    g diff "$prev_commit" "$curr_commit" > "$patchfile"
    g log -1 --pretty=%B "$curr_commit" > "$commitmsg"
    
    prev_commit="$curr_commit"
    i=$((i + 1))
done

gst

if [ -f ~/base_commit.txt ]; then
    git checkout bench_branch
    git reset --hard "$(cat ~/base_commit.txt)"
else
    echo "Missing base_commit.txt"
    exit 1
fi

for i in $(seq "$num_patches" -1 1); do
    lower=$(( i - 1 ))
    patchfile=~/${i}-${lower}.diff
    commitmsg=~/${i}-${lower}.commit
    
    if [ -s "$patchfile" ]; then
        #patch -p1 < "$patchfile" || { echo "Failed to apply $patchfile"; exit 1; }
        g apply "$patchfile" || { echo "Failed to apply $patchfile"; exit 1; }

        gst

        gaa
        gci --author="A U Thor <author@example.com>" -F "$commitmsg" || { echo "Failed to commit with $commitmsg"; exit 1; }
    else
        echo "Patch file $patchfile is empty, skipping commit."
    fi
done

gst
