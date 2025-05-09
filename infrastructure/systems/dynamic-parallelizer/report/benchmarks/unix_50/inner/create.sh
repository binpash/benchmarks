#!/bin/bash

# Base paths
BASE_DIR="/users/gliargko/dynamic-parallelizer/report/benchmarks/unix_50"
FULL_RUN="$BASE_DIR/full/run"

# Check if the full run file exists
if [[ ! -f "$FULL_RUN" ]]; then
    echo "Error: $FULL_RUN not found."
    exit 1
fi

# Copy and modify the run script for each benchmark directory
for BENCHMARK_NO in {1..36}; do
    TARGET_DIR="$BASE_DIR/$BENCHMARK_NO"
    TARGET_RUN="$TARGET_DIR/run"
    
    # Create the target directory if it doesn't exist
    mkdir -p "$TARGET_DIR"

    # Copy the full run file to the target directory
    cp "$FULL_RUN" "$TARGET_RUN"

    # Determine the input file variable based on the benchmark number
    INPUT_FILE="IN${BENCHMARK_NO}=\$IN/1G-${BENCHMARK_NO}.txt"

    # Replace BENCHMARK_NO, SCRIPT_NAME, and IN environment variables in the copied run file
    sed -i "s/BENCHMARK_NO=\"full\"/BENCHMARK_NO=\"$BENCHMARK_NO\"/" "$TARGET_RUN"
    sed -i "s/SCRIPT_NAME=\"full.sh\"/SCRIPT_NAME=\"$BENCHMARK_NO.sh\"/" "$TARGET_RUN"
    sed -i "s|ENV_VARS=(\"IN=\$RESOURCE_DIR\")|ENV_VARS=(\"IN=\$RESOURCE_DIR\" \"$INPUT_FILE\")|" "$TARGET_RUN"

    # Print a message for confirmation
    echo "Created $TARGET_RUN with BENCHMARK_NO=$BENCHMARK_NO and input file $INPUT_FILE"
done