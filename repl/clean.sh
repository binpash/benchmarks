#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done
KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/repl"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"

CHROMIUM_DIR="${input_dir}/chromium"

cd "$CHROMIUM_DIR" || exit 1
git stash
git checkout main
git branch -D bench_branch 2>/dev/null || true

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
fi

rm "${eval_dir}"/vps-audit-negate-report.txt
rm "${eval_dir}"/vps-audit-negate-processed.txt
rm "${eval_dir}"/vps-audit-report.txt
rm "${eval_dir}"/vps-audit-processed.txt
