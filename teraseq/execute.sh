#!/bin/bash --posix

SUITE_DIR="""$(realpath "$(dirname "$0")")"
export SUITE_DIR

export TIMEFORMAT=%R
cd "$SUITE_DIR" || exit 1

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="teraseq"

script_names="run_dRNASeq.sh"

while IFS= read -r script; do
    script_file="./scripts/$script"
    BENCHMARK_SCRIPT="$(realpath "$script_file")"
    export BENCHMARK_SCRIPT

    echo "$script"
    $KOALA_SHELL "$script_file"
    echo "$?"
done <<< "$script_names"
