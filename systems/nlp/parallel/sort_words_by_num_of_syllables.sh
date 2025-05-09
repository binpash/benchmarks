#!/bin/bash 

# Parallelized script for processing input files

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    in_dir=$2
    out_dir=$3

    TEMPDIR=$(mktemp -d)
    cat > "${TEMPDIR}/${input}.words"
    tr -sc '[AEIOUaeiou\012]' ' ' < "${TEMPDIR}/${input}.words" | awk '{print NF}' > "${TEMPDIR}/${input}.syl"
    paste "${TEMPDIR}/${input}.syl" "${TEMPDIR}/${input}.words" | sort -nr | sed 5q > "${out_dir}/${input}.out"
    rm -rf "${TEMPDIR}"
}
export -f pure_func

process_file() {
    input=$1
    in_dir=$2
    out_dir=$3

    cat "$in_dir/$input" |
        tr -c 'A-Za-z' '[\n*]' |
        grep -v "^\s*$" |
        sort -u |
        pure_func "$input" "$in_dir" "$out_dir"
}
export -f process_file

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" process_file {} "${IN}" "${OUT}"

echo "done"
