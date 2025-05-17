#!/usr/bin/env bash
set -x

TOP="$(git rev-parse --show-toplevel)"

SAMPLE_DIR="$TOP/teraseq/inputs"
outdir="$TOP/teraseq/outputs"
benchmark_dir="$TOP/teraseq"

samples="hsa.dRNASeq.HeLa.polyA.1 hsa.dRNASeq.HeLa.polyA.REL5.1 hsa.dRNASeq.HeLa.polyA.PNK.REL5.1"

echo ">>> MAKE DIRECTORY STRUCTURE <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    mkdir -p "$sdir/logfiles"
    mkdir "$sdir/align"
    mkdir "$sdir/db"
done

echo ">>> CHECK FASTQ <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for $i"

    if [ -f "$sdir/fastq/reads.1.fastq.gz" ]; then
        echo "$sdir/fastq/reads.1.fastq.gz is present, continue."
    else
        echo "$sdir/fastq/reads.1.fastq.gz does not exist, trying to download."
        download="$(grep download "$benchmark_dir/README.md" | grep "$i" | cut -d '|' -f 6 | cut -d '(' -f2  | sed 's/)//' | sed 's#https://##' | tr -d '[:space:]')"
        mkdir -p "$sdir/fastq"
        curl "$download" > "$sdir/fastq/reads.1.fastq.gz"
    fi
done

echo ">>> SETUP DATA <<<"

"$benchmark_dir"/data/run.sh
