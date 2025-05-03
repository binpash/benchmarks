#!/bin/bash

# Find all execute.sh files and calculate LOC for each
results=$(find . -name "execute.sh" -exec wc -l {} + | head -n -1)

# Check if any execute.sh files were found
if [ -z "$results" ]; then
  echo "No execute.sh files found in the directory."
  exit 1
fi

# Process the results to extract LOC and file paths
min_file=$(echo "$results" | sort -n | head -n 1)
max_file=$(echo "$results" | sort -n | tail -n 1)

# Extract line counts and paths
min_loc=$(echo "$min_file" | awk '{print $1}')
min_path=$(echo "$min_file" | awk '{print $2}')

max_loc=$(echo "$max_file" | awk '{print $1}')
max_path=$(echo "$max_file" | awk '{print $2}')

# Print results
echo "Minimum LOC: $min_loc, File: $min_path"
echo "Maximum LOC: $max_loc, File: $max_path"
