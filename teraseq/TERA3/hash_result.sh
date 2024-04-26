#!/bin/bash

source ../PARAMS.sh

samples=(
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
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

