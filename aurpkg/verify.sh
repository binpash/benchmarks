#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath $(dirname "$0"))"
mkdir -p hashes/small

[ ! -d "outputs" ] && echo "Directory 'outputs' does not exist" && exit 1

if [[ "$@" == *"--small"* ]]; then
    hash_folder="hashes/small"
else
    hash_folder="hashes"
fi

if [[ "$@" == *"--generate"* ]]; then
    # Directory to iterate over
    directory="outputs"

    # Loop through all PKGBUILD files in the directory and its subdirectories
    find "$directory" -type f -name "PKGBUILD" | while read -r file
    do
        # Extract the package name from the directory path
        package_name=$(basename "$(dirname "$file")")

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Save the hash to a file with the package name
        echo "$hash" > "$hash_folder/$package_name.hash"

        # Print the filename and hash
        echo "$hash_folder/$package_name.hash $hash"
    done

    exit 0
fi

# Loop through all PKGBUILD files in the directory and its subdirectories
find "outputs" -type f -name "PKGBUILD" | while read -r file
do
    # Extract the package name from the directory path
    package_name=$(basename "$(dirname "$file")")

    if [ ! -f "$hash_folder/$package_name.hash" ]; then
        echo "Hash file for $package_name does not exist."
        continue
    fi

    # Generate SHA-256 hash
    hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

    # Read the stored hash
    stored_hash=$(cat "$hash_folder/$package_name.hash")

    if [ "$hash" == "$stored_hash" ]; then
        echo "$package_name: OK"
    else
        echo "$package_name: FAILED"
    fi
done
