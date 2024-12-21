#!/bin/bash

cd "$(realpath $(dirname "$0"))"

HASH_FILE="hashes.txt"

# Create or clear the hash file
> "$HASH_FILE"

echo "Generating hashes for .run files..."

for file in outputs/*.run; do
  if [[ -f "$file" ]]; then
    # Calculate the SHA-256 hash and append to the hash file
    hash=$(shasum -a 256 "$file" | awk '{ print $1 }')
    echo "$(basename "$file") $hash" >> "$HASH_FILE"
    echo "Hash for $(basename "$file") added."
  fi
done

echo "Hashes saved to $HASH_FILE."

