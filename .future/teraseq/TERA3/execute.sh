#!/bin/sh
#
# Run TERA3 preprocessing, alignment, and postprocessing
#

cd /root/TERA-Seq_manuscript/samples # dliu

. ../PARAMS.sh

threads=6
assembly="hg38"

####################################################################################################

samples="hsa.dRNASeq.HeLa.total.REL3.1 hsa.dRNASeq.HeLa.total.REL3.2 hsa.dRNASeq.HeLa.total.REL3.3"

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

echo ">>> REMOVE REL3 ADAPTOR <<<"

# First, we take only max. last 200 to identify reads with adapter and avoid internal trimming
# With REL3 it's more complicated than with REL5 because first we have RTA (and the Guppy RTA trimming doesn't work well)

if [ -z "$CONDA_PREFIX" ]; then
    echo "Variable \$CONDA_PREFIX is not set. Please make sure you specified if in PARAMS.sh."
    exit
fi

. "$CONDA_PREFIX"/bin/activate # Source Conda base
conda activate teraseq

. "$INSTALL"/cutadapt-2.5/venv/bin/activate

# Identify reads with adapter from last 200 bp of each read
N=200 # for rel3
for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    # Trim read up to max. last 200 bp
    gunzip -c "$sdir"/fastq/reads.1.sanitize.fastq.gz | seqkit seq -m $((N+1)) \
    | seqkit subseq -r -$N:-1 | gzip -c > "$sdir"/fastq/tmp1.fastq.gz

    # Get shorter reads just not to exclude them at the start
    gunzip -c "$sdir"/fastq/reads.1.sanitize.fastq.gz | seqkit seq -M $N \
    | gzip -c > "$sdir"/fastq/tmp2.fastq.gz

    cat "$sdir"/fastq/tmp1.fastq.gz "$sdir"/fastq/tmp2.fastq.gz > "$sdir"/fastq/tmp.fastq.gz && rm "$sdir"/fastq/tmp1.fastq.gz "$sdir"/fastq/tmp2.fastq.gz

    cutadapt \
        -a GTGTCAGTCACTTCCA \
        --overlap 16 \
        --minimum-length 0 \
        --error-rate 0.18 \
        --untrimmed-output "$sdir"/fastq/tmp.wo_rel3.fastq.gz \
        --output "$sdir"/fastq/tmp.w_rel3.fastq.gz \
        "$sdir"/fastq/tmp.fastq.gz \
        > "$sdir"/logfiles/cutadapt.rel3.last$N.log 2>&1 # dliu remove &
done
# dliu remove wait

# Get only reads reads with adapter and run trimming once more
for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    zcat "$sdir"/fastq/tmp.w_rel3.fastq.gz | paste - - - - | cut -f1 | sed 's/^@//g' \
        > "$sdir"/fastq/reads.1.sanitize.w_rel3.names.txt # dliu remove &
    zcat "$sdir"/fastq/tmp.wo_rel3.fastq.gz | paste - - - - | cut -f1 | sed 's/^@//g' \
        > "$sdir"/fastq/reads.1.sanitize.wo_rel3.names.txt # dliu remove &
    # dliu remove wait

    rm "$sdir"/fastq/tmp.fastq.gz "$sdir"/fastq/tmp.w_rel3.fastq.gz "$sdir"/fastq/tmp.wo_rel3.fastq.gz

    zcat "$sdir"/fastq/reads.1.sanitize.fastq.gz | seqtk subseq - "$sdir"/fastq/reads.1.sanitize.w_rel3.names.txt \
        | gzip -c > "$sdir"/fastq/tmp.w_rel3.fastq.gz # dliu remove &
    zcat "$sdir"/fastq/reads.1.sanitize.fastq.gz | seqtk subseq - "$sdir"/fastq/reads.1.sanitize.wo_rel3.names.txt \
        | seqkit seq -m 25 | gzip -c > "$sdir"/fastq/reads.1.sanitize.wo_rel3.fastq.gz # dliu remove &

    # dliu remove wait

    cutadapt \
        -a GTGTCAGTCACTTCCA \
        --overlap 16 \
        --minimum-length 25 \
        --error-rate 0.18 \
        --untrimmed-output "$sdir"/fastq/ishouldbeempty.fastq.gz \
        --output "$sdir"/fastq/reads.1.sanitize.w_rel3.fastq.gz \
        "$sdir"/fastq/tmp.w_rel3.fastq.gz \
        > "$sdir"/logfiles/cutadapt.rel3.log 2>&1

    rm "$sdir"/fastq/tmp.w_rel3.fastq.gz
done
# dliu remove wait

deactivate

echo ">>> MERGE READS WITH AND WITHOUT REL3 ADAPTOR <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    cat "$sdir"/fastq/reads.1.sanitize.w_rel3.fastq.gz "$sdir"/fastq/reads.1.sanitize.wo_rel3.fastq.gz \
        > "$sdir"/fastq/reads.1.sanitize.rel3_trim.fastq.gz # dliu remove &
done
# dliu remove wait

echo ">>> ALIGN READS TO RIBOSOMAL (ALL ENSEMBL + SILVA-HUMAN) <<<"

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
        "$DATA_DIR"/$assembly/minimap2.17/ensembl-transcripts-wRibo.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.rel3_trim.fastq.gz \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam

    samtools view -H "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam > \
        "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam
    samtools view -@ $threads -F4 "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam \
        | grep -v -P "\tENST" >> "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam
    rm "$sdir"/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam

    cat "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam | cut -f1 | sort | uniq > \
        "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.reads.txt # dliu remove &

    samtools view -@ $threads -bh "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam | samtools sort -@ $threads - > \
        "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.bam

    # dliu remove wait
    rm "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.sam
done
# dliu remove wait

echo ">>> EXTRACT NON-RIBOSOMAL READS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    seqkit grep -nvf "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.reads.txt \
        "$sdir"/fastq/reads.1.sanitize.rel3_trim.fastq.gz \
        -o "$sdir"/fastq/reads.1.sanitize.noribo.fastq.gz # dliu remove &

    seqkit grep -nvf "$sdir"/align/reads.1.sanitize.toRibosomal.sorted.reads.txt \
        "$sdir"/fastq/reads.1.sanitize.fastq.gz \
        -o "$sdir"/fastq/reads.1.sanitize.noribo-nanopolish.fastq.gz # dliu remove &
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
        "$DATA_DIR"/$assembly/minimap2.17/transcripts.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.noribo.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam
done
# dliu remove wait

echo ">>> ALIGN READS TO TRANSCRIPTOME (WITH SECONDARY) - TOTAL RNA <<<"

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
        "$DATA_DIR"/$assembly/minimap2.17/transcripts-total.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.noribo.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-total.sorted.bam

    ln -s "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-total.sorted.bam \
    "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam
done
# dliu remove wait

echo ">>> ALIGN READS TO TRANSCRIPTOME (WITH SECONDARY) - TOTAL RNA - FOR NANOPOLISH <<<"

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
        "$DATA_DIR"/$assembly/minimap2.17/transcripts-total.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.noribo-nanopolish.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.noribo-nanopolish.toTranscriptome-total.sorted.bam
done
# dliu remove wait

echo ">>> ALIGN READS TO GENOME (WITH SECONDARY) <<<"

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
        "$DATA_DIR"/$assembly/minimap2.17/genome.k12.mmi \
        "$sdir"/fastq/reads.1.sanitize.rel3_trim.fastq.gz \
    | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
    | grep -v "SA:Z:" \
    | sam-count-secondary --tag X0 \
    | samtools view -b - \
    | samtools sort - \
    > "$sdir"/align/reads.1.sanitize.toGenome.sorted.bam
done
# dliu remove wait

echo ">>> INDEX ALIGNMENTS <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    samtools index \
        "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam # dliu remove &
    samtools index \
        "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-total.sorted.bam # dliu remove &
    ln -s "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome-total.sorted.bam.bai \
        "$sdir"/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam.bai
    samtools index \
        "$sdir"/align/reads.1.sanitize.noribo-nanopolish.toTranscriptome-total.sorted.bam # dliu remove &
    samtools index \
        "$sdir"/align/reads.1.sanitize.toGenome.sorted.bam # dliu remove &
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

echo ">>> ANNOTATE WITH GENIC ELEMENTS (TRANSCRIPTOME) <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_file \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --a_file "$DATA_DIR"/$assembly/genic_elements.mrna.bed \
        --column coding_transcript # dliu remove &
done
# dliu remove wait

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    clipseqtools-preprocess annotate_with_file \
        --database "$sdir"/db/sqlite.db \
        --table transcr \
        --a_file "$DATA_DIR"/$assembly/genic_elements.ncrna.bed \
        --column noncoding_transcript # dliu remove &
done
# dliu remove wait

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

deactivate

echo ">>> ANNOTATE WITH REL3 <<<"

conda activate teraseq

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    echo ">> ANNOTATE WITH REL3 (TRANSCRIPTOME) <<"

    annotate-sqlite-with-fastq \
        --database "$sdir"/db/sqlite.db \
        --db_col_bind "qname" \
        --db_col_add "rel3" \
        --db_tables "transcr" \
        --ifile "$sdir"/fastq/reads.1.sanitize.w_rel3.fastq.gz

    echo ">> ANNOTATE W/WO REL3 (GENOME) <<"

    annotate-sqlite-with-fastq \
        --database "$sdir"/db/sqlite.db \
        --db_col_bind "qname" \
        --db_col_add "rel3" \
        --db_tables "genome" \
        --ifile "$sdir"/fastq/reads.1.sanitize.w_rel3.fastq.gz
done

conda deactivate # dliu match conda activate

# Note: If you have the fast5 files you can continue with the analysis. Check https://github.com/mourelatos-lab/TERA-Seq_manuscript/samples/README.md for the location where to download them.

# echo ">>> NANOPOLISH POLYA <<<"
#
# for i in $samples; do
#     sdir=$SAMPLE_DIR/$i
#     echo " Working for" $i
#
#     if [ -d "$sdir/fast5" ]; then
# 	    if [ `find $sdir/fast5 -maxdepth 1 -type f -name '*.fast5' | wc -l` != 0 ]; then
#             nanopolish index \
#                 --directory $sdir/fast5/ \
#                 $sdir/fastq/reads.1.fastq.gz
#
#             nanopolish polya \
#                 --reads $sdir/fastq/reads.1.fastq.gz \
#                 --bam $sdir/align/reads.1.sanitize.noribo-nanopolish.toTranscriptome-total.sorted.bam \
#                 --genome $DATA_DIR/$assembly/transcripts-total.fa \
#                 --threads $threads \
#                 > $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab
#
#             tail -n+2 $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab | egrep -w "PASS|NOREGION" \
#                 | cut -f1,9 > $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab
#         else
#             echo "It seems that $sdir/fast5 directory is empty. Please check you downloaded and uncompressed fast5 tar.gz archive and placed the files in $sdir/fast5."
#         fi
#     else
#         echo "It seems that $sdir/fast5 directory doesn't exist. Please check you created $sdir/fast5, downloaded and uncompressed fast5 tar.gz and placed the files in $sdir/fast5."
#     fi
# done
#
# echo ">>> ANNOTATE WITH POLYA (GENOME) <<<"
#
# for i in $samples; do
#     sdir=$SAMPLE_DIR/$i
#     echo " Working for" $i
#
#     if [ -f "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" ]; then
#         annotate-sqlite-with-file \
#             --db_col_add "polya" --db_col_bind "qname" \
#             --db_tables "genome" --database "$sdir/db/sqlite.db" \
#             --ifile "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" \
#             --round # dliu remove &
#     else
#         echo "It seems $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab doesn't exist. Not going to annotate the db with poly(A) tail lengths."
#     fi
# done
# # dliu remove wait
#
# echo ">>> ANNOTATE WITH POLYA (TRANSCRIPTOME) <<<"
#
# for i in $samples; do
#     sdir=$SAMPLE_DIR/$i
#     echo " Working for" $i
#
#     if [ -f "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" ]; then
#         annotate-sqlite-with-file \
#             --db_col_add "polya" --db_col_bind "qname" \
#             --db_tables "transcr" --database "$sdir/db/sqlite.db" \
#             --ifile "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" \
#             --round # dliu remove &
#     else
#         echo "It seems $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab doesn't exist. Not going to annotate the db with poly(A) tail lengths."
#     fi
# done
# # dliu remove wait

echo ">>> ALL DONE <<<"
