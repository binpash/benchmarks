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
    for d in outputs/*/; do
        b=$(basename "$d")
        out="$hash_folder/$b.hashes"
        : > "$out"
        find "$d" -type f | sort | while read -r f; do
            rel=${f#outputs/}
            printf '%s  %s\n' "$(shasum -a 256 "$f" | awk '{print $1}')" "$rel" >> "$out"
        done
    done
    exit
fi

mismatch=0
for d in outputs/*/; do
    b=$(basename "$d")
    ref="$hash_folder/$b.hashes"
    [[ -f $ref ]] || { echo "$b missing reference"; mismatch=1; continue; }
    tmp=$(mktemp)
    find "$d" -type f | sort | while read -r f; do
        rel=${f#outputs/}
        printf '%s  %s\n' "$(shasum -a 256 "$f" | awk '{print $1}')" "$rel" >> "$tmp"
    done
    diff -u "$ref" "$tmp" | grep -E '^[+-][0-9a-f]' && mismatch=1
    rm -f "$tmp"
done

echo nlp $mismatch