#!/bin/sh
# set -e
# tag: count_consonant_sequences

mkdir -p "$OUTPUT_DIR"

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | tr '[a-z]' '[A-Z]' | tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | sort | uniq -c > $OUTPUT_DIR/$input.out
done
