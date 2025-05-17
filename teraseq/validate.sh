#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

hash_folder="hashes"
output_folder="outputs"

generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
    --min) hash_folder="hashes/min" ;;
    --small) hash_folder="hashes/small" ;;
    esac
done

if $generate; then
    find "$output_folder" -type f | sort | xargs md5sum > "$hash_folder/outputs.hashes"
    exit 0
fi

mismatch=0
tmpfile=$(mktemp)
find "$output_folder" -type f | sort | xargs md5sum > "$tmpfile"
if ! diff -q "$hash_folder/outputs.hashes" "$tmpfile" > /dev/null; then
    diff -u "$hash_folder/outputs.hashes" "$tmpfile"
    mismatch=1
fi

echo "teraseq $mismatch"
