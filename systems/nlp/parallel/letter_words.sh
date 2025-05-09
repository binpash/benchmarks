#!/bin/bash

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

process_file() {
    input=$1
    in_dir=$2
    out_dir=$3

    cat "$in_dir/$input" |
        tr -c 'A-Za-z' '[\n*]' |
        grep -v "^\s*$" |
        grep -c '^....$' > "${out_dir}/${input}.out0"

    cat "$in_dir/$input" |
        tr -c 'A-Za-z' '[\n*]' |
        grep -v "^\s*$" |
        sort -u |
        grep -c '^....$' > "${out_dir}/${input}.out1"
}
export -f process_file

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" process_file {} "${IN}" "${OUT}"

echo "done"