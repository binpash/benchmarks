#!/usr/bin/env bash

cd "$(dirname "$0")" || exit 1

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="web-index"

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

for arg in "$@"; do
    case "$arg" in
        --small) INPUT_FILE="$INPUTS/index_small.txt" ;;
        --min) INPUT_FILE="$INPUTS/index_min.txt" ;;
        *) INPUT_FILE="$INPUTS/index.txt" ;;
    esac
done

mkdir -p "$OUTPUT_BASE"

echo "web-index"
BENCHMARK_SCRIPT="$(realpath "./scripts/ngrams.sh")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$INPUT_FILE")"
export BENCHMARK_INPUT_FILE

"$BENCHMARK_SHELL" "./scripts/ngrams.sh" "$OUTPUT_BASE"
echo $?
