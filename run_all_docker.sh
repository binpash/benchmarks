#!/bin/bash
BENCHMARKS_DIR="/benchmarks"
EXCLUDE_DIR="infrastructure"
SCRIPT_NAME="main.sh"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

args=("$@")

run_benchmarks() {
    for BENCH in "$BENCHMARKS_DIR"/*/; do
        BENCH_NAME=$(basename "$BENCH")

        if [ "$BENCH_NAME" = "$EXCLUDE_DIR" ]; then
            continue
        fi

        echo "Running benchmark: $BENCH_NAME"
        "$BENCHMARKS_DIR/$SCRIPT_NAME" "$BENCH_NAME" "${args[@]}" || echo "Benchmark $BENCH_NAME failed!"
    done
}

echo "Running benchmarks inside Docker container..."
chmod +x "$BENCHMARKS_DIR/$SCRIPT_NAME"
run_benchmarks
