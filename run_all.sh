#!/bin/bash

BENCHMARKS_DIR="/benchmarks"
EXCLUDE_DIR="infrastructure"
SCRIPT_NAME="main.sh"
IMAGE_NAME="koala"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

run_locally=false
args=()

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --bare)
            run_locally=true
            ;;
        *)
            args+=("$arg")
            ;;
    esac
done

run_benchmarks() {
    for BENCH in "$BENCHMARKS_DIR"/*/; do
        BENCH_NAME=$(basename "$BENCH")

        if [[ "$BENCH_NAME" == "$EXCLUDE_DIR" ]]; then
            continue
        fi

        echo "Running benchmark: $BENCH_NAME"
        "$BENCHMARKS_DIR/$SCRIPT_NAME" "$BENCH_NAME" "${args[@]}" || echo "Benchmark $BENCH_NAME failed!"
    done
}

if $run_locally; then
    echo "Running benchmarks locally..."
    chmod +x "$BENCHMARKS_DIR/$SCRIPT_NAME"
    run_benchmarks
else
    echo "Running benchmarks inside Docker..."
    
    docker build -t "$IMAGE_NAME" .

    docker run --cap-add NET_ADMIN --cap-add NET_RAW -it "$IMAGE_NAME" "$BENCHMARK_SHELL" -c "
        chmod +x /benchmarks/$SCRIPT_NAME
        $(declare -f run_benchmarks)  # Export function to the container
        run_benchmarks \"${args[@]}\"
    "
fi
