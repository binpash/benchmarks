#!/bin/bash 

# Parallelized script for generating trigrams

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/4_3b/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    in_dir=$2
    out_dir=$3

    TEMPDIR=$(mktemp -d)

    # Process the input
    cat "$in_dir/$input" |
        tr -c 'A-Za-z' '[\n*]' |
        grep -v "^\s*$" > "${TEMPDIR}/${input}.words"

    tail +2 "${TEMPDIR}/${input}.words" > "${TEMPDIR}/${input}.nextwords"
    tail +2 "${TEMPDIR}/${input}.words" > "${TEMPDIR}/${input}.nextwords2"

    paste "${TEMPDIR}/${input}.words" "${TEMPDIR}/${input}.nextwords" "${TEMPDIR}/${input}.nextwords2" |
        sort |
        uniq -c > "${out_dir}/${input}.trigrams"

    # Cleanup
    rm -rf "${TEMPDIR}"
}
export -f pure_func

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" pure_func {} "${IN}" "${OUT}"

echo "done"
