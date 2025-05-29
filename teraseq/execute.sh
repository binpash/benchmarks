#!/bin/bash --posix

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=full ;; # small uses a subset of full inputs
    --min) size=min ;;
    esac
done

export SIZE="$size" # for PARAMS.sh

SUITE_DIR="""$(realpath "$(dirname "$0")")"
export SUITE_DIR

export TIMEFORMAT=%R
cd "$SUITE_DIR" || exit 1

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="teraseq"

script_names="data
run_dRNASeq
run_5TERA"

if [[ "$size" == "small" ]]; then
    script_names="run_5TERA-short"
fi

while IFS= read -r script; do
    script_file="./scripts/$script.sh"
    BENCHMARK_SCRIPT="$(realpath "$script_file")"
    export BENCHMARK_SCRIPT

    echo "$script"
    $KOALA_SHELL "$script_file"
    echo "$?"
done <<< "$script_names"
