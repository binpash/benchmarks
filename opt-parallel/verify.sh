#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath "$(dirname "$0")")" || exit 1

[ ! -d "outputs" ] && echo "Directory 'outputs' does not exist" && exit 1

for arg in "$@"; do
    case "$arg" in
    --small) hash_folder="hashes/small" ;;
    --min) hash_folder="hashes/min" ;;
    esac
done

mkdir -p "$hash_folder"

if [[ " $* " == *" --generate "* ]]; then
    for file in outputs/*.out; do
        filename=$(basename "$file" .out)
        hash=$(shasum -a 256 "$file" | awk '{print $1}')
        echo "$hash" > "$hash_folder/$filename.hash"
        echo "$hash_folder/$filename.hash $hash"
    done
    exit 0
fi

for file in outputs/*.out; do
    filename=$(basename "$file" .out)
    hash=$(shasum -a 256 "$file" | awk '{ print $1 }')
    echo "$hash" > "outputs/$filename.hash"
    diff "$hash_folder/$filename.hash" "outputs/$filename.hash" > /dev/null
    match="$?"
    echo "$filename $match"
done
