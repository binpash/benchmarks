#!/bin/bash

set -e
set -o pipefail

REPO_TOP=$(git rev-parse --show-toplevel)
EXCLUDE_DIR="infrastructure"
SCRIPT_NAME="main.sh"
KOALA_SHELL=${KOALA_SHELL:-bash}
if [[ "$1" =~ ^- ]]; then
    OUTPUT_PATH="${REPO_TOP}/infrastructure/dynamic_analysis"
else
    OUTPUT_PATH="$1"
    shift
fi

args=("$@" --resources)

pushd "$REPO_TOP"/infrastructure

sudo apt-get update && apt-get install -y \
    git \
    python3 python3-pip python3-venv \
    build-essential libtool m4 automake cloc autoconf time gawk
python3 -m venv venv && . venv/bin/activate && pip install -r requirements.txt


checked_command() {
    local description="$1"
    shift
    echo "Generating $description" >&2
    if ! "$@"; then
	echo "Error generating $description, static analysis failed" >&2
	exit 1
    fi
}

echo "Creating target directory (if it doesn't exist)"
mkdir -p target plots

checked_command "target/scripts_to_benchmark.csv" python3 scripts_to_benchmark.py | sort > target/scripts_to_benchmark.csv

checked_command "target/lines_of_code.csv" python3 count_lines_of_code.py | sort > target/lines_of_code.csv

checked_command "target/nodes_in_scripts.csv" python3 count_nodes_in_scripts.py | sort > target/nodes_in_scripts.csv

echo "Generating target/shellmetrics.sh"
if ! wget --quiet -O target/shellmetrics.sh https://raw.githubusercontent.com/shellspec/shellmetrics/b3bfff2af6880443112cdbf2ea449440b30ab9b0/shellmetrics; then
  echo "Error downloading target/shellmetrics.sh"
  exit 1
fi

checked_command "target/shellmetrics.sh" chmod a+x target/shellmetrics.sh

checked_command "target/cyclomatic.csv" python3 get_cyclomatic.py | sort > target/cyclomatic.csv

checked_command "plots/koala-stx-analysis.pdf" python3 viz/syntax.py plots

echo "Static analysis complete, all targets generated successfully. Syntax analysis heatmap is located at infrastructure/plots/koala-stx-analysis.pdf"
popd


mkdir -p "$OUTPUT_PATH"

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
rm -f "$REPO_TOP"/infrastructure/target/*.csv

pushd "$REPO_TOP"/infrastructure

checked_command_dynamic() {
    local description="$1"
    shift
    echo "Generating $description" >&2
    if ! "$@"; then
	echo "Error generating $description, dynamic analysis failed" >&2
	exit 1
    fi
}


checked_command_dynamic "target/dynamic_analysis.jsonl" python3 dynamic_analysis.py | sort > target/dynamic_analysis.jsonl

python3 viz/dynamic.py "$OUTPUT_PATH"
cat "$OUTPUT_PATH/benchmark_stats.txt"
echo "Dynamic analysis complete, target generated successfully. Dynamic analysis plots are located in $OUTPUT_PATH"

popd