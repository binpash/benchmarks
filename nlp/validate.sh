#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath "$(dirname "$0")")" || exit 1

hash_folder="hashes"

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

mkdir -p "$hash_folder"

if $generate; then
    for dir in outputs/*/; do
        script=$(basename "$dir")
        out="$hash_folder/$script.hashes"
        : > "$out"
        find "$dir" -type f | grep -v '\.hash$' | sort | while read -r f; do
            rel=${f#outputs/}
            printf '%s  %s\n' "$(shasum -a 256 "$f" | awk '{print $1}')" "$rel" >> "$out"
        done
    done
    exit
fi

mismatch=0
for dir in outputs/*/; do
    script=$(basename "$dir")
    ref="$hash_folder/$script.hashes"
    [[ -f $ref ]] || { echo "$script missing reference"; mismatch=1; continue; }
    tmp=$(mktemp)
    find "$dir" -type f | grep -v '\.hash$' | sort | while read -r f; do
        rel=${f#outputs/}
        printf '%s  %s\n' "$(shasum -a 256 "$f" | awk '{print $1}')" "$rel" >> "$tmp"
    done
    if ! diff -q "$ref" "$tmp" > /dev/null; then
        echo "Mismatch in $script:"
        diff -u "$ref" "$tmp"
        mismatch=1
    fi
    rm -f "$tmp"
done

echo "nlp $mismatch"
