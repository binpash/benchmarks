#!/bin/bash

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.2_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

process_file() {
    input=$1
    in_dir=$2
    out_dir=$3

    cat "$in_dir/$input" |
        tr -c 'A-Za-z' '[\n*]' |
        grep -v "^\s*$" |
        tr -sc 'AEIOUaeiou' '[\012*]' |
        sort |
        uniq -c |
        awk '$1 >= 1000' > "${out_dir}/${input}.out"
}
export -f process_file

# Use GNU Parallel to process files concurrently
ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" process_file {} "${IN}" "${OUT}"

echo "done"