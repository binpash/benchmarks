#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/web-index/outputs"
hashes_dir="${REPO_TOP}/web-index/hashes"

# create hashes directory if it does not exist
if [ ! -d "${hashes_dir}" ]; then
    mkdir "${hashes_dir}"
fi

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
    # generate hashes and store in hashes directory for all *grams.txt files
    for file in $(find ${eval_dir} -name "*grams.txt"); do
        hash=$(md5sum <(sort "$file") | cut -d ' ' -f 1)
        echo "${hash}" > "${hashes_dir}/$(basename ${file})${suffix}.hash"
        echo "$file $hash"
    done
    exit 0
fi

# verify hashes for all *grams.txt files
for file in $(find ${eval_dir} -name "*grams.txt"); do
    hash=$(md5sum <(sort "$file") | cut -d ' ' -f 1)
    expected_hash=$(cat "${hashes_dir}/$(basename ${file})${suffix}.hash")
    match=0
    if [[ "${hash}" != "${expected_hash}" ]]; then
        match=1
    fi
    echo "$file $match"
done