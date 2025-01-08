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
    for file in "$directory"/*.out
    do
        # Extract the filename without the directory path and extension
        filename=$(basename "$file" .out)

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Save the hash to a file
        echo "$hash" > "$hash_folder/$filename.hash"

        # Print the filename and hash
        echo "$hash_folder/$filename.hash $hash"
    done

    exit 0
fi

# Loop through all directories in the parent directory
for file in "outputs"/*.out
do
    # Extract the filename without the directory path and extension
    filename=$(basename "$file" .out)

    # Generate SHA-256 hash
    shasum -a 256 "$file" | awk '{ print $1 }' > "$file.hash"

    # Compare the hash with the hash in the hashes directory
    diff "$hash_folder/$filename.hash" "$file.hash" > /dev/null
    match="$?"

    # Print the filename and hash
    echo "outputs/$filename.out $match"
done
