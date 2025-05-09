#!/bin/sh
# tag: four-letter words
# set -e

mkdir -p "$OUTPUT_DIR"

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -c '^....$' > $OUTPUT_DIR/$input.out0
    cat $INPUT_FILE/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | grep -c '^....$'  > $OUTPUT_DIR/$input.out1
done
