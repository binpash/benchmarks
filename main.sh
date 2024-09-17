#!/usr/bin/env bash

# Main benchmark driver.
# This will be responsible for:
# 0. Starting up a containerized environment
# 1. Preparing benchamrk inputs
# 2. Running the benchmark
# 3. Collecting the results
# 4. Cleaning up the environment
# 5. Reporting the results

# Set the Docker image name
IMAGE_NAME="shell-benchmark"
BENCHMARK_TOP=$(git rev-parse --show-toplevel)

cd "$BENCHMARK_TOP" || exit 1

# Function to run the benchmark by passing the correct benchmark folder
run_benchmark() {
    local benchmark_name=$1
    local benchmark_folder="$PWD/$benchmark_name"

    echo "Running benchmark: $benchmark_name"

    docker run --rm -e BENCHMARK_FOLDER="$benchmark_folder" -v "$PWD/$benchmark_name:/benchmarks/$benchmark_name" "$IMAGE_NAME"
    
    echo "Benchmark $benchmark_name completed."
}

# Check if a benchmark was passed as an argument
if [ -z "$1" ]; then
    echo "No benchmarks specified. Please provide benchmark names as arguments."
    echo "Example: ./driver.sh benchmark1 benchmark2"
    exit 1
fi

# Loop through all the benchmarks passed as arguments
for benchmark in "$@"; do
    run_benchmark "$benchmark"
done

echo "All benchmarks completed."
