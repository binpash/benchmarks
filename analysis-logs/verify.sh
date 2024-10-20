#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-logs"
input_dir="${eval_dir}/input"
hashes_dir="${eval_dir}/hashes"
results_dir="${eval_dir}/results"
mkdir -p $results_dir

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd
if [[ "$@" == *"--generate"* ]]; then
    md5sum results$suffix/* > $hashes_dir/results$suffix.md5sum
fi

okay=0
if ! md5sum --check --quiet $hashes_dir/results$suffix.md5sum; then
    okay=1
    echo "img_convert $suffix failed verification"
fi
exit $okay
