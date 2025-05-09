#!/bin/sh
# verses with 2 or more, 3 or more, exactly 2 instances of light.

mkdir -p "$OUTPUT_DIR"

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | grep -c 'light.\*light'                                 > $OUTPUT_DIR/$input.out0
    cat $INPUT_FILE/$input | grep -c 'light.\*light.\*light'                         > $OUTPUT_DIR/$input.out1
    cat $INPUT_FILE/$input | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > $OUTPUT_DIR/$input.out2
done
