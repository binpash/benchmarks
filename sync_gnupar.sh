#!/bin/bash

ORIGINAL_REPO="/benchmarks"
MODIFIED_SUITE="./temp/benchmarks/infrastructure/systems/GNU-Parallel"

sync_benchmarks() {
    local src="$1"
    local dest="$2"

    for benchmark in "$src"/*; do
        if [ -d "$benchmark" ]; then
            benchmark_name=$(basename "$benchmark")

            if [ "$benchmark_name" == "infrastructure" ]; then
                echo "Skipping infrastructure directory"
                continue
            fi

            dest_benchmark="$dest/$benchmark_name"

            mkdir -p "$dest_benchmark"

            rsync -av --exclude 'input' "$benchmark/" "$dest_benchmark/"

            echo "Synced $benchmark_name to $dest_benchmark"
        fi
    done
}

sync_benchmarks "$MODIFIED_SUITE" "$ORIGINAL_REPO"