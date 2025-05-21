#!/bin/bash --posix

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=full ;; # small uses a subset of full inputs
    --min) size=min ;;
    esac
done

TOP=$(realpath "$(dirname "$0")")

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="prog-inf"

SUITE_DIR="""$(realpath "$(dirname "$0")")"
export SUITE_DIR

export TIMEFORMAT=%R
cd "$SUITE_DIR" || exit 1

mkdir -p "outputs"
export INDEX="$TOP/inputs/index.$size.txt"

script_file="$TOP/scripts/proginf.sh"
BENCHMARK_SCRIPT=$(realpath "$script_file")
export BENCHMARK_SCRIPT
echo "proginf.sh"
$KOALA_SHELL "$script_file"
echo "$?"
