#!/bin/bash
# tag: verse_2om_3om_2instances
# set -e
# verses with 2 or more, 3 or more, exactly 2 instances of light.

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_7/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    # Process the file
    grep -c 'light.*light' < "$IN/$input" > "${OUT}/${input}.out0" &
    grep -c 'light.*light.*light' < "$IN/$input" > "${OUT}/${input}.out1" &
    grep 'light.*light' < "$IN/$input" | grep -vc 'light.*light.*light' > "${OUT}/${input}.out2" &
done

wait

echo 'done';
# rm -rf ${OUT}
