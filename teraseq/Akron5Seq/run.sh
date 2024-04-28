#!/bin/sh
#
# Run Akron5-Seq (SRR6360508) preprocessing, alignment, and postprocessing
#

cd /root/TERA-Seq_manuscript/samples # dliu

. ../PARAMS.sh

threads=6
assembly="hg38"

####################################################################################################

samples="hsa.Akron5Seq.HeLa.whole.2"

echo ">>> REMOVE ALL ADAPTORS (PAIRED-END) <<<"

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
        -a GTGTCAGTCACTTCCAGCGG \
        -A CCGCATCGTCCTCCCT \
        -u 16 \
        --overlap 5 \
        --minimum-length 25 \
        --error-rate 0.15 \
        --output "$sdir"/fastq/reads.1.noadapt.fastq.gz \
        --paired-output "$sdir"/fastq/reads.2.noadapt.fastq.gz \
        "$sdir"/fastq/reads.1.fastq.gz \
        "$sdir"/fastq/reads.2.fastq.gz \
        > "$sdir"/logfiles/cutadapt.log 2>&1 # dliu remove &
done
# dliu remove wait

deactivate

echo ">>> ALIGN READS (PAIRED-END) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    STAR \
        --genomeDir "$DATA_DIR"/$assembly/STAR-2.7.2b/ \
        --readFilesIn "$sdir"/fastq/reads.1.noadapt.fastq.gz "$sdir"/fastq/reads.2.noadapt.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN $threads \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All \
        --outFilterMultimapScoreRange 0 \
        --alignIntronMax 50000 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outFilterMatchNmin 8 \
        --outFilterMatchNminOverLread 0.6 \
        --outFileNamePrefix "$sdir"/align/reads. \
        --sjdbOverhang 100 \
        --alignSJDBoverhangMin 1 \
        --sjdbGTFfile "$DATA_DIR"/$assembly/genes.gtf \
        --quantMode TranscriptomeSAM \
        --outReadsUnmapped Fastx \
        > "$sdir"/logfiles/star.log 2>&1
done

echo ">>> CLEANUP AND GZIP FILES <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    gzip -f "$sdir"/align/reads.Unmapped.out.mate1 # dliu remove &
    gzip -f "$sdir"/align/reads.Unmapped.out.mate2 # dliu remove &
done
# dliu remove wait

CONDA_PATH=$CONDA_PREFIX # Temporary store path to the Conda environment
conda deactivate

echo ">>> MARK DUPLICATES <<<"

. "$INSTALL"/perl-virtualenv/teraseq/bin/activate

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    "$CONDA_PATH"/bin/samtools view -h \
        "$sdir"/align/reads.Aligned.sortedByCoord.out.bam \
    | sam-mark-dups --primary-tag "XD:Z" --paired \
    | "$CONDA_PATH"/bin/samtools view -b - \
    > "$sdir"/align/reads.Aligned.sortedByCoord.marked.out.bam # dliu remove &
done
# dliu remove wait

deactivate

echo ">>> SAM TO FASTQ FOR PROPER NON DUPLICATE PAIRS <<<"

conda activate teraseq

PICARD_RUN="java -jar -Xmx8g $(ls "$CONDA_PREFIX"/share/picard-*/picard.jar)"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools view -h -f 2 -F 1024 \
        "$sdir"/align/reads.Aligned.sortedByCoord.marked.out.bam \
    | $PICARD_RUN SamToFastq \
        INPUT=/dev/stdin \
        FASTQ="$sdir"/fastq/reads.1.noadapt.aligned.nodups.fastq \
        SECOND_END_FASTQ="$sdir"/fastq/reads.2.noadapt.aligned.nodups.fastq \
        QUIET=true

    gzip -f "$sdir"/fastq/reads.1.noadapt.aligned.nodups.fastq # dliu remove &
    gzip -f "$sdir"/fastq/reads.2.noadapt.aligned.nodups.fastq # dliu remove &
done
# dliu remove wait

echo ">>> ALIGN NON DUPLICATE READ 1 <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    STAR \
        --genomeDir "$DATA_DIR"/$assembly/STAR-2.7.2b/ \
        --readFilesIn "$sdir"/fastq/reads.1.noadapt.aligned.nodups.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN $threads \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All \
        --outFilterMultimapScoreRange 0 \
        --alignIntronMax 50000 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outFilterMatchNmin 8 \
        --outFilterMatchNminOverLread 0.7 \
        --outFileNamePrefix "$sdir"/align/reads.1. \
        --sjdbOverhang 80 \
        --alignSJDBoverhangMin 1 \
        --seedSearchStartLmax 15 \
        --sjdbGTFfile "$DATA_DIR"/$assembly/genes.gtf \
        --quantMode TranscriptomeSAM \
        > "$sdir"/logfiles/star.nodups.log 2>&1
done
# dliu remove wait

echo ">>> SORT TRANSCRIPTOME ALIGNMENT <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools sort --threads $threads \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.out.bam \
    > "$sdir"/align/reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam
    rm "$sdir"/align/reads.1.Aligned.toTranscriptome.out.bam
done
# dliu remove wait

echo ">>> HOMOGENIZE FILENAMES <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    ln -sf \
        reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.final.bam
    ln -sf \
        reads.1.Aligned.sortedByCoord.out.bam \
        "$sdir"/align/reads.1.Aligned.final.bam
done
# dliu remove wait

echo ">>> INDEX ALIGNMENTS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools index \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.final.bam
    samtools index \
        "$sdir"/align/reads.1.Aligned.final.bam
done
# dliu remove wait

conda deactivate

echo ">>> SAM TO SQLITE (TRANSCRIPTOME) <<<"

. "$INSTALL"/perl-virtualenv/teraseq/bin/activate

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/align/reads.1.Aligned.toTranscriptome.final.bam \
    | "$CONDA_PATH"/bin/samtools view -h -F 272 - \
    | sam_to_sqlite-short \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --drop # dliu remove &
done
# dliu remove wait

echo ">>> ANNOTATE WITH GENIC ELEMENTS (TRANSCRIPTOME) <<<"

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

echo ">>> ANNOTATE WITH GENIC ELEMENTS (GENOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/align/reads.1.Aligned.final.bam \
    | "$CONDA_PATH"/bin/samtools view -h -F 256 - \
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
