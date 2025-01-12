#!/usr/bin/env bash

cd "$(dirname "$0")"

export BENCHMARK_CATEGORY="web-index"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
directory_path="inputs/articles"

if [ ! -d "$directory_path" ]; then
	    echo "Error: Directory does not exist."
      exit 1
fi

# ensure a local ./tmp directory exists for sorting
mkdir -p ./tmp
export TMPDIR=$PWD/tmp

INPUTS="$PWD/inputs"
OUTPUT_BASE="$PWD/outputs/ngrams"

if [[ "$@" == *"--small"* ]]; then
    export INPUT_FILE="$INPUTS/index_small.txt"
else
    export INPUT_FILE="$INPUTS/index.txt"
fi

mkdir -p "$OUTPUT_BASE"

echo "web-index"
export BENCHMARK_SCRIPT="$(realpath "./scripts/ngrams.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$INPUT_FILE")"
time $BENCHMARK_SHELL ./scripts/ngrams.sh "$OUTPUT_BASE"
echo $?
