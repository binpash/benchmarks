#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath $(dirname "$0"))"

hash_folder="hashes"

mkdir -p $hash_folder

if [[ "$@" == *"--generate"* ]]; then
    directory="output"

    # get total number of files
    total_files=$(ls "$directory"/*.bam | wc -l)
    current_file=0

    # Loop through all .bam files in the directory
    for file in "$directory"/*.bam
    do
        current_file=$((current_file + 1))
        echo "Processing file $current_file of $total_files..."

        # Extract the filename without the directory path and extension
        filename=$(basename "$file" .bam)

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Save the hash to a file
        echo "$hash" > "$hash_folder/$filename.hash"

        # Print the filename and hash
        echo "File: $hash_folder/$filename.hash | SHA-256 Hash: $hash"
    done
fi

# Loop through all directories in the parent directory
for folder in "output"
do
    # Remove trailing slash
    folder=${folder%/}

    echo "Verifying folder: $folder"

    # total number of files
    total_files=$(ls "$folder"/*.bam | wc -l)
    current_file=0

    # Loop through all .bam files in the current directory
    for file in "$folder"/*.bam
    do
        current_file=$((current_file + 1))
        # echo "Processing file $current_file of $total_files..."

        # Extract the filename without the directory path and extension
        filename=$(basename "$file" .bam)

        if [ ! -f "$folder/$filename.hash" ]; then
            # Generate SHA-256 hash
            hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

            # Save the hash to a file
            echo "$hash" > "$folder/$filename.hash"
        fi

        # Compare the hash with the hash in the hashes directory
        diff "$hash_folder/$filename.hash" "$folder/$filename.hash" > /dev/null
        match="$?"

        # Print the filename and match
        echo "$folder/$filename $match"
    done
done
exit 0
