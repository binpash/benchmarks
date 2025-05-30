#!/bin/bash
TOP=$(git rev-parse --show-toplevel)
EXCLUDE_DIR="${TOP}/infrastructure"
KOALA_SHELL=${KOALA_SHELL:-bash}
if [[ "$1" =~ ^- ]]; then
    OUTPUT_PATH=/tmp/plots/dynamic_analysis
else
    OUTPUT_PATH="$1"
    shift
fi
mkdir -p "$OUTPUT_PATH"

args=("$@" --resources)

run_benchmarks() {
    for BENCH in "$TOP"/*/; do
        BENCH_NAME=$(basename "$BENCH")

        if [ "$BENCH_NAME" = "$EXCLUDE_DIR" ]; then
            continue
        fi

        echo "Running benchmark: $BENCH_NAME"
        $TOP/main.sh "$BENCH_NAME" "${args[@]}" || echo "Benchmark $BENCH_NAME failed!"
    done
}

echo "Running benchmarks inside Docker container..."
chmod +x "$TOP/$SCRIPT_NAME"
run_benchmarks

rm -f "$TOP"/infrastructure/target/dynamic_analysis.jsonl
rm -f "$TOP"/infrastructure/target/*.csv
cd "$TOP/infrastructure" || exit 1
make

python3 "$TOP/infrastructure/viz/dynamic.py" "$OUTPUT_PATH"
cat "$OUTPUT_PATH/benchmark_stats.txt"
echo '--------------------------------------------'
echo "Dynamic analysis plots saved to $OUTPUT_PATH"
echo '--------------------------------------------'
