#!/bin/bash

source ../PARAMS.sh

samples=(
    "hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1"
)

for f in "${samples[@]}"; do
    sdir="$SAMPLE_DIR/$f"
    for g in find "$sdir"; do
        echo -n "$g: "
        if [ -f "$g" ]
        then
            if [[ "$g" == *.gz ]]
            then
                zcat "$g" | md5sum
            else
                md5sum < "$g"
            fi
        fi
    done
done

