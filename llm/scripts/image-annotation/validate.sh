#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/llm/scripts/image-annotation"
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"
source "$REPO_TOP/venv/bin/activate"

suffix=""
generate=false
for arg in "$@"; do
    case "$arg" in
        --generate) generate=true ;;
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done
hashes_dir="${hashes_dir}/jpg$suffix"
mkdir -p "$hashes_dir"

outputs_dir="$outputs_dir/jpg$suffix"
mkdir -p "$outputs_dir"

if $generate; then
    cd $outputs_dir || exit 1
    bench=image-annotation$suffix
    md5sum $outputs_dir/* > "$hashes_dir/$bench.md5sum"
    exit 0
fi

cd $outputs_dir || exit 1
bench=image-annotation$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?