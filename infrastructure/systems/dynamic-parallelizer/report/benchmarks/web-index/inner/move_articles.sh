#!/bin/bash

# Define the base directory
base_directory="input1000"

# Check if the base directory exists
if [ ! -d "$base_directory" ]; then
    echo "Base directory does not exist: $base_directory"
    exit 1
fi

# Navigate to the base directory
cd "$base_directory"

# Create a tar archive of the en/articles directory
tar -czvf en_articles.tar.gz en/articles

echo "Archive created: $(pwd)/en_articles.tar.gz"
