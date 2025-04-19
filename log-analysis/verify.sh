#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
hashes_dir="${eval_dir}/hashes"
outputs_dir="${eval_dir}/outputs"
mkdir -p "${outputs_dir}"

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

cd "$outputs_dir" || exit # md5sum computes paths relative to cd

if $generate; then
    md5sum "pcaps$suffix"/* > "$hashes_dir/pcaps$suffix.md5sum"
    md5sum "nginx$suffix"/* > "$hashes_dir/nginx$suffix.md5sum"
    exit 0
fi

bench=pcaps$suffix
md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
echo $bench $?

bench=nginx$suffix
md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
echo $bench $?
