#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath $(dirname "$0"))"

mkdir -p hashes/small

if [[ "$@" == *"--small"* ]]; then
    hash_folder="hashes/small"
else
    hash_folder="hashes"
fi

if [[ "$@" == *"--generate"* ]]; then
    # Directory to iterate over
    directory="outputs"

    # Loop through all .out files in the directory
    find "$directory" -type f -name '*.out' | while read -r file;
    do
        # Extract the filename and dirname
        filename=$(basename "$file" .out)
        dirname=$(dirname "${file#$directory/}")

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Create subdirectory if not already
        mkdir -p $hash_folder/$dirname

        # Save the hash to a file
        echo "$hash" > "$hash_folder/$dirname/$filename.hash"

        # Print the filename and hash
        echo "$hash_folder/$dirname/$filename.hash" "$hash"
    done

    exit 0
fi

# Loop through all directories in the parent directory
for folder in "outputs"/*
do
    # Loop through all .out files in the current directory
    find "$folder" -type f -name '*.out' | while read -r file;
    do
        # Extract the filename and dirname
        filename=$(basename "$file" .out)
        dirname=$(basename "$(dirname "$file")") # is the script_name

        # Generate SHA-256 hash
        shasum -a 256 "$file" | awk '{ print $1 }' > "$file.hash"

        # Compare the hash with the hash in the hashes directory
        diff "$hash_folder/$dirname/$filename.hash" "$file.hash" > /dev/null
        match="$?"

        # Print the filename and match
        echo "$dirname/$filename $match"
    done
done
