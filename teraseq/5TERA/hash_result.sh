#!/bin/bash

source ../PARAMS.sh

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
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

