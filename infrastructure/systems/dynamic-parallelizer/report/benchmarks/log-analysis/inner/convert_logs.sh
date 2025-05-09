#!/bin/bash

# Define the directory containing the log files and mappings
export PATH=$PATH:$HOME/.local/bin
export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

download_dir="$PASH_SPEC_TOP/report/resources/log-analysis"
benchmark_dir="$PASH_SPEC_TOP/report/benchmarks/log-analysis"
log_dir="$download_dir"
mapping_file="${benchmark_dir}/object_mappings.sort"
tool="${benchmark_dir}/recreate"

# Ensure the mappings file exists
if [ ! -f "$mapping_file" ]; then
    echo "Error: Mappings file $mapping_file not found!"
    exit 1
fi

# Ensure the tool exists
if [ ! -x "$tool" ]; then
    echo "Error: Tool $tool not found or not executable!"
    exit 1
fi

# Process each uncompressed log file in the directory
for file in "$log_dir"/*; do
    if [ -f "$file" ] && [[ "$file" == *.gz ]]; then
        output_file="${file%.gz}.log"
        echo "Processing $file and saving to $output_file"
    	gzip -dc "$file" | "$tool" "$mapping_file" > "$output_file"
    fi
done
