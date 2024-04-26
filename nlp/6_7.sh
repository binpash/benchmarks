#!/bin/bash
# tag: verse_2om_3om_2instances
# set -e
# verses with 2 or more, 3 or more, exactly 2 instances of light.

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/6_7}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | grep -c 'light.\*light'                                 > ${OUT}/${input}.out0
    cat $IN/$input | grep -c 'light.\*light.\*light'                         > ${OUT}/${input}.out1
    cat $IN/$input | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > ${OUT}/${input}.out2
done

echo 'done';
