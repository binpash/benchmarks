#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/input"
scripts_dir="${eval_dir}/scripts"
hashes_dir="${eval_dir}/hashes"
results_dir="${eval_dir}/results"
mkdir -p $results_dir

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd

if [[ "$@" == *"--generate"* ]]; then
    md5sum pcaps$suffix/* > $hashes_dir/pcaps$suffix.md5sum
    md5sum nginx$suffix/* > $hashes_dir/nginx$suffix.md5sum
    exit 0
fi

bench=pcaps$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum 
echo $bench $?

bench=nginx$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?
