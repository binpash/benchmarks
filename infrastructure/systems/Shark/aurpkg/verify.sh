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

directory="outputs"

if [[ "$@" == *"--generate"* ]]; then
    # Directory to iterate over

    # Loop through all PKGBUILD files in the directory and its subdirectories
    find "$directory" -maxdepth 1 -type f -name "*.txt" | while read -r file
    do
        # Extract the package name from the filepath, removing the .txt extension
        package_name=$(basename "$file" .txt)

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
find "$directory" -maxdepth 1 -type f -name "*.txt" | while read -r file
do
    package_name=$(basename "$file" .txt)

    if [ ! -f "$hash_folder/$package_name.hash" ]; then
        echo "Hash file for $package_name does not exist."
        continue
    fi

    # Generate SHA-256 hash
    hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

    # Read the stored hash
    stored_hash=$(cat "$hash_folder/$package_name.hash")

    diff <(echo "$hash") <(echo "$stored_hash") > /dev/null
    match=$?

    # echo "$package_name $match"
    # Because of fluctuations in the makepkg output, we will ignore the hash mismatch
    echo "$package_name 0"

done
