#!/bin/bash

# Parallelized script for processing input files

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    in_dir=$2
    out_dir=$3
    category=$4

    TEMPDIR=$(mktemp -d)
    tr -sc '[A-Z][a-z]' '[\012*]' > "${TEMPDIR}/${input}.words"
    tail +2 "${TEMPDIR}/${input}.words" > "${TEMPDIR}/${input}.nextwords"
    tail +3 "${TEMPDIR}/${input}.words" > "${TEMPDIR}/${input}.nextwords2"
    paste "${TEMPDIR}/${input}.words" "${TEMPDIR}/${input}.nextwords" "${TEMPDIR}/${input}.nextwords2" | sort | uniq -c | sort -nr | sed 5q > "${out_dir}/${input}.${category}.out"
    rm -rf "${TEMPDIR}"
}
export -f pure_func

process_file() {
    input=$1
    in_dir=$2
    out_dir=$3

    cat "$in_dir/$input" | grep 'the land of' | pure_func "$input" "$in_dir" "$out_dir" "0"
    cat "$in_dir/$input" | grep 'And he said' | pure_func "$input" "$in_dir" "$out_dir" "1"
}
export -f process_file

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" process_file {} "${IN}" "${OUT}"

echo "done"