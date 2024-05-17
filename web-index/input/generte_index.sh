#!/bin/bash
if [ $# -eq 0 ]; then
	    echo "Usage: $0 <directory_path>"
	        exit 1
fi

# Directory path is the first argument
directory_path=$1

# Check if the directory exists
if [ ! -d "$directory_path" ]; then
	    echo "Error: Directory does not exist."
	        exit 1
fi

# Ensure a local ./tmp directory exists for sorting
mkdir -p ./tmp
export TMPDIR=./tmp

# Find all files, remove prefix, sort them, and write to a text file
find "$directory_path" -type f | sed 's|./wikipedia/en/articles/||' | sort > index.txt

echo "File paths have been saved to all_files_paths.txt"

