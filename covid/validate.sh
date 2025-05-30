#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"

suffix=""
generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
    --small) suffix="_small" ;;
    --min) suffix="_min" ;;
    esac
done

if $generate; then
    mkdir -p "$hashes_dir"
    # give relative paths to md5sum
    (
        cd "$outputs_dir" || exit 1
        md5sum "outputs$suffix"/* >"$hashes_dir/outputs$suffix.md5sum"
    )
    exit 0
fi

# give relative paths to md5sum
(
    cd "$outputs_dir" || exit 1
    md5sum --check --quiet --status "$hashes_dir/outputs$suffix.md5sum"
)
echo covid-mts$suffix $?
