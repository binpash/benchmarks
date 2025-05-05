#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath "$(dirname "$0")")" || exit 1
mkdir -p hashes/small

if [ ! -d "outputs" ]; then
    echo "Directory 'outputs' does not exist"
    exit 1
fi

hash_folder="hashes"
generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
        --small) hash_folder="hashes/small" ;;
        --min) hash_folder="hashes/min" ;;
    esac
done

mkdir -p "$hash_folder"

if $generate; then
    # Directory to iterate over
    directory="outputs"

    # Loop through all .out files in the directory
    for file in "$directory"/*.out; do
        # Extract the filename without the directory path and extension
        filename=$(basename "$file" .out)

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Save the hash to a file
        echo "$hash" >"$hash_folder/$filename.hash"

        # Print the filename and hash
        echo "$hash_folder/$filename.hash $hash"
    done

    exit 0
fi

# Loop through all directories in the parent directory
for file in outputs/*.out; do
    # Extract the filename without the directory path and extension
    filename=$(basename "$file" .out)

    # Generate SHA-256 hash
    hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

    # Save the hash to a file
    echo "$hash" > "outputs/$filename.hash"

    diff "$hash_folder/$filename.hash" "outputs/$filename.hash" > /dev/null
    match="$?"
    # Print the filename and hash
    echo "outputs/$filename $match"
done
