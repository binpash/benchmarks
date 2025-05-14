#!/bin/sh
#
# Run Mouse SIRV (PRJEB27590) preprocessing, alignment, and postprocessing
#

cd /root/TERA-Seq_manuscript/samples # dliu

. ../PARAMS.sh

threads=6
assembly="mm10"

####################################################################################################

samples="mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680375.1 mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680379.1 mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363657.1 mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363659.1"

echo ">>> SANITIZE FASTQ HEADERS <<<"

. "$INSTALL"/perl-virtualenv/teraseq/bin/activate

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    zcat "$sdir"/fastq/reads.1.fastq.gz \
    | fastq-sanitize-header --input - --delim : --keep 0 \
    | gzip \
    > "$sdir"/fastq/reads.1.sanitize.fastq.gz # dliu remove &
done
# dliu remove wait

deactivate

if [ -z "$CONDA_PREFIX" ]; then
    echo "Variable \$CONDA_PREFIX is not set. Please make sure you specified if in PARAMS.sh."
    exit
fi

. "$CONDA_PREFIX"/bin/activate # Source Conda base
conda activate teraseq

echo ">>> ALIGN READS TO RIBOSOMAL (ALL ENSEMBL + SILVA-MOUSE) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    minimap2 \
        -a \
        -x map-ont \
        -k 12 \
        -p 1 \
        -u f \
        -t $threads \
        --secondary=yes \
        "$DATA_DIR"/$assembly/minimap2.17/ensembl-transcripts-wRibo_sirv1.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.fastq.gz \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam

    samtools view -H "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam > \
        "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam
    samtools view -@ $threads -F4 "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam \
        | grep -v -P "\tENSMUST" | grep -v -P "\tSIRV" >> "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam

    cat "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam | cut -f1 | sort | uniq > \
        "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.reads.txt # dliu remove &

    samtools view -@ $threads -bh "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam | samtools sort -@ $threads - > \
        "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.bam && rm "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam
done
# dliu remove wait

echo ">>> EXTRACT NON-RIBOSOMAL READS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    seqkit grep -nvf "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.reads.txt "$sdir"/fastq/reads.1.sanitize.fastq.gz \
        -o "$sdir"/fastq/reads.1.sanitize.noribo.fastq.gz # dliu remove &
done
# dliu remove wait

echo ">>> ALIGN READS TO TRANSCRIPTOME (WITH SECONDARY) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    minimap2 \
        -a \
        -x map-ont \
        -k 12 \
        -p 1 \
        -u f \
        -t $threads \
        --secondary=yes \
        "$DATA_DIR"/$assembly/minimap2.17/transcripts_sirv1.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.noribo.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam

    ln -s "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam \
        "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam
done
# dliu remove wait

echo ">>> ALIGN READS TO GENOME (WITH SECONDARY) <<<"
# --splice-flank=no doesn't require specific bases to be present around the splice sites and we don't know how it looks here; has to be used in SIRV spike-in and maybe here as well (https://lh3.github.io/minimap2/minimap2.html)

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    minimap2 \
        -a \
        -x splice \
        -k 12 \
        -p 1 \
        -u b \
        -t $threads \
        --secondary=yes \
        "$DATA_DIR"/$assembly/minimap2.17/genome_sirv1.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.toGenome.sorted.bam
done
# dliu remove wait

echo ">>> ALIGN READS TO GENOME (WITH SECONDARY) - FOR SIRV <<<"
# --splice-flank=no doesn't require specific bases to be present around the splice sites and we don't know how it looks here; has to be used in SIRV spike-in and maybe here as well (https://lh3.github.io/minimap2/minimap2.html)

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    minimap2 \
        -a \
        -x splice \
        -k 12 \
        -p 1 \
        -u b \
        -t $threads \
        --secondary=yes \
        --splice-flank=no \
        "$DATA_DIR"/$assembly/minimap2.17/genome_sirv1.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.toGenome.forSIRV.sorted.bam
done
# dliu remove wait

echo ">>> INDEX ALIGNMENTS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools index \
        "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam # dliu remove &
    ln -s "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam.bai \
        "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam.bai
    samtools index \
        "$sdir"/align/reads.1.sanitize.toGenome.sorted.bam # dliu remove &
    samtools index \
        "$sdir"/align/reads.1.sanitize.toGenome.forSIRV.sorted.bam # dliu remove &
    # dliu remove wait
done

CONDA_PATH=$CONDA_PREFIX # Temporary store path to the Conda environment
conda deactivate

echo ">>> SAM TO SQLITE (TRANSCRIPTOME) <<<"

. "$INSTALL"/perl-virtualenv/teraseq/bin/activate

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
    | "$CONDA_PATH"/bin/samtools view -h -F 4 -F 16 -F 2048 - \
    | sam_to_sqlite \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --records_class GenOOx::Data::File::SAMminimap2::Record \
        --drop # dliu remove &
done
# dliu remove wait

echo ">>> SAM TO SQLITE (GENOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/align/reads.1.sanitize.toGenome.sorted.bam \
    | "$CONDA_PATH"/bin/samtools view -h -F 4 -F 2048 - \
    | sam_to_sqlite \
        --database "$sdir"/db/sqlite.db \
        --table genome \
        --records_class GenOOx::Data::File::SAMminimap2::Record \
        --drop # dliu remove &
done
# dliu remove wait

echo ">>> ANNOTATE WITH GENIC ELEMENTS (GENOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_genic_elements \
        --database "$sdir"/db/sqlite.db \
        --table genome \
        --gtf "$DATA_DIR"/$assembly/genes-polya.gtf # dliu remove &
done
# dliu remove wait

echo ">>> ALL DONE <<<"
