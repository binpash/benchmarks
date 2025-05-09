#!/bin/bash
# Parallelized script for processing with two input sources

IN=${IN:-$SUITE_DIR/inputs/pg}
INPUT2=${INPUT2:-$SUITE_DIR/inputs/exodus}
OUT=${1:-$SUITE_DIR/outputs/8.3_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    input_dir=$2
    input2=$3
    out_dir=$4
    TEMPDIR=$(mktemp -d)

    # Process input1
    cat "$input_dir/$input" | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u > "${TEMPDIR}/${input}1.types"

    # Process input2
    cat "$input2" | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > "${TEMPDIR}/${input}2.types"

    # Combine and count types
    sort "${TEMPDIR}/${input}1.types" "${TEMPDIR}/${input}2.types" "${TEMPDIR}/${input}2.types" | uniq -c | head > "${out_dir}/${input}.out"

    # Cleanup
    rm -rf "${TEMPDIR}"
}
export -f pure_func

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" pure_func {} "${IN}" "${INPUT2}" "${OUT}"

echo "done"
