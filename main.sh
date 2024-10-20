BENCHMARK="$1"

cd "$(dirname "$0")/$BENCHMARK" || exit 1

error()
{
    echo "Error: $1" > /dev/stderr
    exit 1
}

# Download dependencies
./deps.sh || error "Failed to download dependencies for $BENCHMARK"

# Fetch inputs
./input.sh || error "Failed to fetch inputs for $BENCHMARK"

# Run benchmark
./run.sh > $BENCHMARK.out 2> $BENCHMARK.err || error "Failed to run $BENCHMARK"

# Verify output
./verify.sh > $BENCHMARK.hash || error "Failed to verify output for $BENCHMARK"
