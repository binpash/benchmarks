#!/bin/sh
#
# Run Ribo-Seq (SRR3306589) preprocessing, alignment, and postprocessing
#

cd /root/TERA-Seq_manuscript/samples # dliu

. ../PARAMS.sh

threads=6
assembly="hg38"

####################################################################################################

samples="hsa.RiboSeq.HeLa.async.2"

echo ">>> TRIM 3' ADAPTOR <<<"
if [ -z "$CONDA_PREFIX" ]; then
    echo "Variable \$CONDA_PREFIX is not set. Please make sure you specified if in PARAMS.sh."
    exit
fi

. "$CONDA_PREFIX"/bin/activate # Source Conda base
conda activate teraseq

. "$INSTALL"/cutadapt-2.5/venv/bin/activate

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cutadapt \
        -a TGGAATTCTCGGGTGCCAAGG \
        --overlap 6 \
        --minimum-length 15 \
        --error-rate 0.1 \
        --discard-untrimmed \
        -o "$sdir"/fastq/reads.1.adtrim.fastq.gz \
        "$sdir"/fastq/reads.1.fastq.gz \
        > "$sdir"/logfiles/cutadapt.log 2>&1 # dliu remove &
done
# dliu remove wait

deactivate

echo ">>> ALIGN READS (WITH STAR) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    STAR \
        --genomeDir "$DATA_DIR"/$assembly/STAR-2.7.2b/ \
        --readFilesIn "$sdir"/fastq/reads.1.adtrim.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN $threads \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All \
        --outFilterMultimapScoreRange 0 \
        --alignIntronMax 50000 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outFilterMatchNmin 15 \
        --outFilterMatchNminOverLread 0.9 \
        --outFileNamePrefix "$sdir"/align/reads.1. \
        --sjdbOverhang 100 \
        --sjdbGTFfile "$DATA_DIR"/$assembly/genes.gtf \
        --quantMode TranscriptomeSAM \
        > "$sdir"/logfiles/star.log 2>&1
done

echo ">>> SORT TRANSCRIPTOME ALIGNMENT <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools sort -@ $threads \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.out.bam \
    > "$sdir"/align/reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam
    rm "$sdir"/align/reads.1.Aligned.toTranscriptome.out.bam
done

echo ">>> CLEAN rRNAs and tRNAs <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    bedtools intersect -s -v \
        -abam "$sdir"/align/reads.1.Aligned.sortedByCoord.out.bam \
        -b "$DATA_DIR"/$assembly/rRNA_tRNA.bed \
    > "$sdir"/align/reads.1.Aligned.sortedByCoord.clean.out.bam # dliu remove &
done
# dliu remove wait

CONDA_PATH=$CONDA_PREFIX # Temporary store path to the Conda environment
conda deactivate

echo ">>> MARK DUPLICATES <<<"

. "$INSTALL"/perl-virtualenv/teraseq/bin/activate

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    "$CONDA_PATH"/bin/samtools view -@ $threads -h \
        "$sdir"/align/reads.1.Aligned.sortedByCoord.clean.out.bam \
    | sam-mark-dups --primary-tag "XD:Z" \
    | "$CONDA_PATH"/bin/samtools view -@ $threads -b - \
    > "$sdir"/align/reads.1.Aligned.sortedByCoord.clean.marked.out.bam
done

echo ">>> MARK DUPLICATES (TRANSCRIPTOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    "$CONDA_PATH"/bin/samtools view -@ $threads -h \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam \
    | sam-mark-dups --primary-tag "XD:Z" \
    | "$CONDA_PATH"/bin/samtools view -@ $threads -b - \
    > "$sdir"/align/reads.1.Aligned.toTranscriptome.sortedByCoord.marked.out.bam
done

echo ">>> HOMOGENIZE FILENAMES <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    ln -sf \
        reads.1.Aligned.toTranscriptome.sortedByCoord.marked.out.bam \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.final.bam
    ln -sf \
        reads.1.Aligned.sortedByCoord.clean.marked.out.bam \
        "$sdir"/align/reads.1.Aligned.final.bam
done

echo ">>> INDEX ALIGNMENTS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    "$CONDA_PATH"/bin/samtools index \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.final.bam
    "$CONDA_PATH"/bin/samtools index \
        "$sdir"/align/reads.1.Aligned.final.bam
done
# dliu remove wait

echo ">>> SAM TO SQLITE (TRANSCRIPTOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/align/reads.1.Aligned.toTranscriptome.final.bam \
    | "$CONDA_PATH"/bin/samtools view -@ $threads -h -F 272 - \
    | sam-mark-dups-add-count-tag \
        --primary-tag "XD:Z" \
        --count-tag "XC:i" \
    | "$CONDA_PATH"/bin/samtools view -@ $threads -h -F 1024 - \
    | sam_to_sqlite-short \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --drop # dliu remove &
done
# dliu remove wait

echo ">>> ANNOTATE WITH GENIC ELEMENTS) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_file \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --a_file "$DATA_DIR"/$assembly/genic_elements.utr5.bed \
        --column utr5 # dliu remove &
done
# dliu remove wait

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_file \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --a_file "$DATA_DIR"/$assembly/genic_elements.cds.bed \
        --column cds # dliu remove &
done
# dliu remove wait

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_file \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --a_file "$DATA_DIR"/$assembly/genic_elements.utr3.bed \
        --column utr3 # dliu remove &
done
# dliu remove wait

echo ">>> SAM TO SQLITE (GENOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/align/reads.1.Aligned.final.bam \
    | "$CONDA_PATH"/bin/samtools view -@ $threads -h -F 256 - \
    | sam-mark-dups-add-count-tag \
        --primary-tag "XD:Z" \
        --count-tag "XC:i" \
    | "$CONDA_PATH"/bin/samtools view -@ $threads -h -F 1024 - \
    | sam_to_sqlite-short \
        --database "$sdir"/db/sqlite.db \
        --table genome \
        --drop # dliu remove &
done
# dliu remove wait

echo ">>> ANNOTATE WITH GENIC ELEMENTS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_genic_elements \
        --database "$sdir"/db/sqlite.db \
        --table genome \
        --gtf "$DATA_DIR"/$assembly/genes.gtf # dliu remove &
done
# dliu remove wait

echo ">>> ALL DONE <<<"
