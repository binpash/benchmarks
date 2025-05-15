#!/bin/bash

set -euo pipefail

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/port-scan"
hash_folder="${eval_dir}/hashes"

generate=false
size=full

for arg in "$@"; do
    case "$arg" in
        --small) size="small" ;;
        --min) size="min" ;;
        --generate) generate=true ;;
    esac
done

outputs_dir="$eval_dir/outputs/$size"
hash_folder="$eval_dir/hashes/$size"
mkdir -p "$hash_folder"

as_popularity_file="$outputs_dir/as_popularity.csv"

if $generate; then
    filename=$(basename "$as_popularity_file")
    hash=$(shasum -a 256 "$as_popularity_file" | awk '{ print $1 }')
    echo "$hash" > "$hash_folder/$filename.hash"
    echo "$hash_folder/$filename.hash $hash"
    exit 0
fi

all_ok=0
filename=$(basename "$as_popularity_file")
hash_file="$hash_folder/$filename.hash"

if [[ ! -f "$as_popularity_file" ]]; then
    echo "Missing data file: $as_popularity_file"
    exit 1
fi

actual_hash=$(shasum -a 256 "$as_popularity_file" | awk '{ print $1 }')

if [[ ! -f "$hash_file" ]]; then
    echo "Missing hash file: $hash_file"
    all_ok=1
else
    expected_hash=$(cat "$hash_file")
    if [[ "$actual_hash" != "$expected_hash" ]]; then
        echo "$filename 1"
        all_ok=1
    else
        echo "$filename 0"
    fi
fi

exit $all_ok
