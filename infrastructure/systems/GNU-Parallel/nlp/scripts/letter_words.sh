#!/bin/bash
# tag: four-letter words
# set -e

# the original script has both versions
# IN=${IN:-$SUITE_DIR/inputs/pg}
# OUT=${1:-$SUITE_DIR/outputs/6_2/}
# ENTRIES=${ENTRIES:-1000}
# mkdir -p "$OUT"

# for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
# do
#     cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -c '^....$' > ${OUT}/${input}.out0
#     cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | grep -c '^....$'  > ${OUT}/${input}.out1
# done

# echo 'done';
# rm -rf "$OUT"


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

rm -f file_list.txt

echo "done"