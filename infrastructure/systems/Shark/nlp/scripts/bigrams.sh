#!/bin/bash 
# tag: bigrams.sh
# set -e

# Bigrams (contrary to our version, this uses intermediary files)
# IN=${IN:-$SUITE_DIR/inputs/pg}
# OUT=${1:-$SUITE_DIR/outputs/4_3/}
# ENTRIES=${ENTRIES:-1000}
# mkdir -p "$OUT"

# pure_func() {
#     input=$1
#     TEMPDIR=$(mktemp -d)
#     cat > ${TEMPDIR}/${input}.input.words
#     tail +2 ${TEMPDIR}/${input}.input.words > ${TEMPDIR}/${input}.input.nextwords
#     paste ${TEMPDIR}/${input}.input.words ${TEMPDIR}/${input}.input.nextwords
#     rm -rf ${TEMPDIR}
# }
# export -f pure_func

# for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
# do
#     cat $IN/$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$"| pure_func $input| sort | uniq -c > ${OUT}/${input}.input.bigrams.out
# done

# echo 'done';
# rm -rf ${OUT}
#!/bin/bash

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/4_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1

    # Generate bigrams directly in memory
    tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | \
    paste - <(tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | tail -n +2)
}

export -f pure_func

# Process each input file
for input in $(find "$IN" -type f | head -n ${ENTRIES}); do
    pure_func "$(basename "$input")" | \
    sort | uniq -c > "${OUT}/$(basename "$input").input.bigrams.out" &
done

# Wait for all background processes to finish
wait

echo 'done'
