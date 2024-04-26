#!/bin/bash

source ../PARAMS.sh

samples=(
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680375.1"
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680379.1"
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363657.1"
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363659.1"
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

