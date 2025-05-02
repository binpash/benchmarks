#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath "$(dirname "$0")")" || exit 1

hash_folder="hashes"
directory="outputs"

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
    # get total number of files
    total_files=$(find "$directory" -maxdepth 1 -name "*.bam" | wc -l)

    # Loop through all .bam files in the directory
    for file in "$directory"/*.bam; do
        # Extract the filename without the directory path and extension
        filename=$(basename "$file" .bam)

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Save the hash to a file
        echo "$hash" > "$hash_folder/$filename.hash"

        # Print the filename and hash
        echo "File: $hash_folder/$filename.hash | SHA-256 Hash: $hash"
    done

    exit 0
fi

# Loop through all .bam files in the current directory
for file in "$directory"/*.bam; do
    # Extract the filename without the directory path and extension
    filename=$(basename "$file" .bam)

    if [ ! -f "$hash_folder/$filename.hash" ]; then
        echo "Error: Hash file for $filename does not exist in $hash_folder."
        echo "Please generate the hash files first using --generate option."
    fi

    # Compare the hash with the hash in the hashes directory
    current_hash=$(shasum -a 256 "$file" | awk '{ print $1 }')
    stored_hash=$(cat "$hash_folder/$filename.hash")

    if [ "$current_hash" = "$stored_hash" ]; then
        match=0
    else
        match=1
    fi

    # Print the filename and match
    echo "$hash_folder/$filename $match"
done
