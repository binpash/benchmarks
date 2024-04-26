#!/bin/bash

source ../PARAMS.sh

samples=(
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.4"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.5"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
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

