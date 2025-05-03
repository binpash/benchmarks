#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
hashes_dir="${eval_dir}/hashes"
outputs_dir="${eval_dir}/outputs"
mkdir -p "${outputs_dir}"

size="full"
generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
        --small) size="small" ;;
        --min) size="min" ;;
    esac
done

cd "$outputs_dir" || exit # md5sum computes paths relative to cd

if $generate; then
    md5sum "pcaps_$size"/* > "$hashes_dir/pcaps_$size.md5sum"
    md5sum "nginx_$size"/* > "$hashes_dir/nginx_$size.md5sum"
    exit 0
fi

bench=pcaps_$size
md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
echo $bench $?

bench=nginx_$size
md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
echo $bench $?
