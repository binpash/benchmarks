#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
EXCLUDE_DIR="infrastructure"
SCRIPT_NAME="main.sh"
KOALA_SHELL=${KOALA_SHELL:-bash}
if [[ "$1" =~ ^- ]]; then
    OUTPUT_PATH="${REPO_TOP}/dynamic_analysis"
else
    OUTPUT_PATH="$1"
    shift
fi
mkdir -p "$OUTPUT_PATH"

args=("$@" --resources)

run_benchmarks() {
    for BENCH in "$REPO_TOP"/*/; do
        BENCH_NAME=$(basename "$BENCH")

        if [ "$BENCH_NAME" = "$EXCLUDE_DIR" ]; then
            continue
        fi

        echo "Running benchmark: $BENCH_NAME"
        $REPO_TOP/$SCRIPT_NAME "$BENCH_NAME" "${args[@]}" || echo "Benchmark $BENCH_NAME failed!"
    done
}

echo "Running benchmarks inside Docker container..."
chmod +x "$REPO_TOP/$SCRIPT_NAME"
run_benchmarks

rm -f "$REPO_TOP"/infrastructure/target/dynamic_analysis.jsonl
cd "$REPO_TOP/infrastructure" || exit 1
make target/dynamic_analysis.jsonl
python3 viz/dynamic.py "$OUTPUT_PATH"
cat "$OUTPUT_PATH/benchmark_stats.txt"
echo "Dynamic analysis plots saved to $OUTPUT_PATH"
