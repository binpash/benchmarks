#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

[ ! -d "outputs" ] && echo "Directory 'outputs' does not exist" && exit 1

size="full"
generate=false
for arg in "$@"; do
    case "$arg" in
        --generate) generate=true ;;
        --small)    size="small"  ;;
        --min)      size="min"    ;;
    esac
done

hash_folder="hashes/$size"
mkdir -p "$hash_folder"
if $generate; then
    find "$output_folder" -type f | sort | xargs md5sum > "$hash_folder/outputs.hashes"
    exit 0
fi
directory="outputs/aurpkg.$size"
input="inputs/packages.$size"

missing=0

while IFS= read -r pkg || [ -n "$pkg" ]; do
    file="$directory/$pkg.txt"
    if [ ! -f "$file" ]; then
        missing=$((missing + 1))
        continue
    fi
    if ! grep -q "Finished making" "$file"; then
        missing=$((missing + 1))
    fi
done < "$input"

if [ "$missing" -eq 0 ]; then
    echo "aurpkg 0"
else
    echo "aurpkg 1"
fi

output_folder="outputs"

mismatch=0
tmpfile=$(mktemp)
find "$output_folder" -type f ! -path "$output_folder/aurpkg.$size/*" | sort | xargs md5sum > "$tmpfile"
if ! diff -q "$hash_folder/outputs.hashes" "$tmpfile" > /dev/null; then
    diff -u "$hash_folder/outputs.hashes" "$tmpfile"
    mismatch=1
fi

echo "prog-inf $mismatch"
