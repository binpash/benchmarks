#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/image_annotation"
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"

suffix=".full"
generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done

if $generate; then
    cd $outputs_dir || exit 1
    bench=image_annotation$suffix
    md5sum $bench/* > "$hashes_dir/$bench.md5sum"
    exit 0
fi

cd $outputs_dir || exit 1
bench=image_annotation$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?