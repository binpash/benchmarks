# !/bin/bash

# Parallelized script for processing light patterns in input files

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_7/}
ENTRIES=${ENTRIES:-1000}

mkdir -p "$OUT"

process_file() {
    input=$1
    in_dir=$2
    out_dir=$3

    cat "$in_dir/$input" | grep -c 'light.\*light'                                 > "${out_dir}/${input}.out0"
    cat "$in_dir/$input" | grep -c 'light.\*light.\*light'                         > "${out_dir}/${input}.out1"
    cat "$in_dir/$input" | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > "${out_dir}/${input}.out2"
}

export -f process_file

# Use GNU Parallel to process files concurrently
ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" process_file {} "${IN}" "${OUT}"

echo 'done'
