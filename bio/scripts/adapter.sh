#!/bin/bash
#
# Run visualization of Cutadapt adapter trimming
#

TOP=$(git rev-parse --show-toplevel)

outdir="$TOP/bio/outputs/teraseq"
scripts="$TOP/bio/scripts"
mkdir -p "$outdir"/fastq

. "$scripts/PARAMS.sh"

adapter_orig="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" # original, full-length REL5.long
ad_name="REL5.long" # adapter name

threads=8
assembly="hg38"

RES_DIR="$outdir/adapter/results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> TRIM ADAPTER - TRANSCRIPTS <<<"

max=${#adapter_orig} # get max length of the adapter if we need to

input="$DATA_DIR/$assembly/transcripts.fa"

sdir="$RES_DIR/transcripts/${ad_name}"
mkdir -p $sdir/logfiles

#for err in $(seq 15 35); do
err=15
while [ "$err" -le 35 ]; do
    min=10
    while [ "$min" -le "$max" ]; do
        rnd=$RANDOM
        # Subset adapter to desired length from the 5p end
        #len=$[$max-$min] # We need to get part of the adapter from the right, not left like in REL
        len=$((max - min))
        adapter=${adapter_orig:$len:$max}

        echo "cutadapt \
            -g $adapter \
            --overlap $min \
            --minimum-length 1 \
            --error-rate 0.$err \
            --output /dev/null \
            $input \
            &> $sdir/logfiles/cutadapt.l${min}.e${err}.log"
        min=$((min + 1))
    done
    err=$((err + 1))
done > $RES_DIR/cmds.txt

cat $RES_DIR/cmds.txt | parallel -j $threads --load 95% --noswap '{}'
rm $RES_DIR/cmds.txt

echo ">>> TRIM ADAPTER - LIBRARIES <<<"

samples="hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"

max=${#adapter_orig} # get max length of the adapter if we need to

adapter=$adapter_orig

adapter=X${adapter} # Add anchoring to the 5' end

for sample in $samples; do
    sdir=$RES_DIR/$sample/cutadapt
    mkdir -p $sdir/logfiles

    err=15
    while [ "$err" -le 35 ]; do
        min=10
        while [ "$min" -le "$max" ]; do
            # The original adapter
            rnd=$RANDOM

            echo "cutadapt \
                -g $adapter \
                --overlap $min \
                --minimum-length 1 \
                --error-rate 0.$err \
                --output /dev/null \
                $SAMPLE_DIR/$sample/fastq/reads.1.fastq.gz \
                &> $sdir/logfiles/cutadapt.orig.l${min}.e${err}.log"

            # Shuffled adapter  - three rounds
            for round in $(seq 1 3); do
                adapter_shuf=$(echo $adapter_orig | fold -w1 | shuf | tr -d '\n') # Shuffle adapter to get random control
                adapter_shuf=X${adapter_shuf}
                rnd=$RANDOM

                echo "cutadapt \
                    -g $adapter_shuf \
                    --overlap $min \
                    --minimum-length 1 \
                    --error-rate 0.$err \
                    --output /dev/null \
                    $SAMPLE_DIR/$sample/fastq/reads.1.fastq.gz \
                    &> $sdir/logfiles/cutadapt.shuf.l${min}.e${err}.${round}.log"
                round=$((round + 1))
            done

            min=$((min + 1))
        done
        err=$((err + 1))
    done
done > $RES_DIR/cmds.txt

cat $RES_DIR/cmds.txt | parallel -j $threads --load 95% --noswap '{}'
rm $RES_DIR/cmds.txt

echo ">>> VISUALIZE TRIMMING - TRANSCRIPTS <<<"

cutadapt-transcripts.R "$RES_DIR/transcripts/${ad_name}/logfiles" \
    "$RES_DIR/transcripts/${ad_name}/cutadapt-heatmap.png"

echo ">>> VISUALIZE TRIMMING - LIBRARIES <<<"

for i in $samples; do
    echo "Working for $i"
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    cutadapt-libraries.R "$sdir/cutadapt/logfiles" \
        "$sdir/cutadapt/cutadapt-heatmap.png" \
        "$RES_DIR/transcripts/${ad_name}/cutadapt-heatmap.tsv"
done

echo ">>> VISUALIZE TRIMMED LENGTHS <<<"

trimming.R "$DATA_DIR/adapter/trimming.tsv" "$RES_DIR/trimming.pdf"

echo ">>> TEST CUTADAPT SENSITIVITY <<<"

mkdir tmp
rmdir tmp

echo ">>> ALL DONE <<<"
