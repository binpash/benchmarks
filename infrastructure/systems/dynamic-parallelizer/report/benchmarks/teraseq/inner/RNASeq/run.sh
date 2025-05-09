#!/bin/sh
#
# Run RNA-Seq (ENCSR000CPR.1) preprocessing, alignment, and postprocessing
#

cd /root/TERA-Seq_manuscript/samples # dliu

. ../PARAMS.sh

threads=6
assembly="hg38"

####################################################################################################

samples="hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1"

echo ">>> TRIM 3' ADAPTOR AND LOW QUALITY ENDS <<<"

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
        -a AGATCGGAAGAGCGGTTCAG \
        -A AGATCGGAAGAGCGTCGTGT \
        --overlap 3 \
        -q 5,5 \
        --trim-n \
        --minimum-length 15 \
        --error-rate 0.1 \
        -o "$sdir"/fastq/reads.1.adtrim.fastq.gz \
        -p "$sdir"/fastq/reads.2.adtrim.fastq.gz \
        "$sdir"/fastq/reads.1.fastq.gz \
        "$sdir"/fastq/reads.2.fastq.gz \
        > "$sdir"/logfiles/cutadapt.log 2>&1 # dliu remove &
done
# dliu remove wait

deactivate

echo ">>> ALIGN READS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    STAR \
        --genomeDir "$DATA_DIR"/$assembly/STAR-2.7.2b-annot/ \
        --readFilesIn "$sdir"/fastq/reads.1.adtrim.fastq.gz "$sdir"/fastq/reads.2.adtrim.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN $threads \
        --outSAMunmapped Within \
        --outSAMattributes All \
        --outSAMheaderHD @HD VN:1.4 SO:coordinate \
        --outSAMattrRGline ID:"${name}" PL:Illumina PU:"${name}" SM:"${name}" \
        --outSAMtype BAM SortedByCoordinate \
        --outMultimapperOrder Random \
        --outFilterMultimapScoreRange 1 \
        --alignMatesGapMax 1000000 \
        --alignIntronMax 1000000 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --alignIntronMin 20 \
        --alignSJoverhangMin 5 \
        --alignSJDBoverhangMin 3 \
        --twopassMode Basic \
        --outFilterType BySJout \
        --outFilterMatchNmin 10 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 1.0 \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterMatchNminOverLread 0.66 \
        --sjdbOverhang 100 \
        --sjdbGTFfile "$DATA_DIR"/$assembly/genes.gtf \
        --quantMode GeneCounts TranscriptomeSAM \
        --outFileNamePrefix "$sdir"/align/reads.12. \
        > "$sdir"/logfiles/star.log 2>&1
done

echo ">>>> SORT TRANSCRIPTOME ALIGNMENT <<<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools sort --threads $threads -T "$sdir"/align/deleteme \
        "$sdir"/align/reads.12.Aligned.toTranscriptome.out.bam \
    > "$sdir"/align/reads.12.Aligned.toTranscriptome.sortedByCoord.out.bam && mv "$sdir"/align/reads.12.Aligned.toTranscriptome.sortedByCoord.out.bam "$sdir"/align/reads.12.Aligned.toTranscriptome.out.bam
done

echo ">>>> HOMOGENIZE ALIGNMENT FILES <<<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    ln -s "$sdir"/align/reads.12.Aligned.sortedByCoord.out.bam \
        "$sdir"/align/reads.1.Aligned.sortedByCoord.out.bam
    ln -s "$sdir"/align/reads.12.Aligned.toTranscriptome.out.bam \
        "$sdir"/align/reads.1.Aligned.toTranscriptome.out.bam
done

echo ">>>> INDEX BAM FILES <<<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools index -@ $threads "$sdir"/align/reads.1.Aligned.sortedByCoord.out.bam
    samtools index -@ $threads "$sdir"/align/reads.1.Aligned.toTranscriptome.out.bam
done

echo ">>> ALL DONE<<<"
