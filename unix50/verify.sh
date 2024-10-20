#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/unix50"
input_dir="${eval_dir}/input"
scripts_dir="${eval_dir}/scripts"
hashes_dir="${eval_dir}/hashes"
results_dir="${eval_dir}/results"

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd

if [[ "$@" == *"--generate"* ]]; then
    md5sum result$suffix/* > $hashes_dir/result$suffix.md5sum
fi

okay=0
if ! md5sum --check --quiet $hashes_dir/result$suffix.md5sum; then
    okay=1
    echo "unix50 $suffix failed verification"
fi

exit $okay
