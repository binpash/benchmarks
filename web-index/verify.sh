#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/web-index"
hashes_dir="${eval_dir}/hashes"

# create hashes directory if it does not exist
if [ ! -d "${hashes_dir}" ]; then
    mkdir "${hashes_dir}"
fi

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

if [[ "$@" == *"--generate"* ]]; then
    # generate hashes and store in hashes directory for all *grams.txt files
    for file in $(find ${eval_dir} -name "*grams.txt"); do
        echo "Generating hash for ${file}"
        hash=$(md5sum ${file} | cut -d ' ' -f 1)
        echo "Hash: ${hash}"
        echo "${hash}" > "${hashes_dir}/$(basename ${file})${suffix}.hash"
    done
    exit 0
fi

# verify hashes for all *grams.txt files
for file in $(find ${eval_dir} -name "*grams.txt"); do
    hash=$(md5sum ${file} | cut -d ' ' -f 1)
    expected_hash=$(cat "${hashes_dir}/$(basename ${file})${suffix}.hash")
    if [[ "${hash}" != "${expected_hash}" ]]; then
        exit 1
    fi
done
exit 0