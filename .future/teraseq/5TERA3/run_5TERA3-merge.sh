#!/bin/bash
#
# Merge all 5TERA3 together to one sample for future analysis
#

trap 'echo Exiting "$BASH_COMMAND" with status $?' EXIT

set -e # dliu

cd /root/TERA-Seq_manuscript/samples # dliu

. ../PARAMS.sh

threads=6

####################################################################################################
if [ -z "$CONDA_PREFIX" ]; then
    echo "Variable \$CONDA_PREFIX is not set. Please make sure you specified if in PARAMS.sh."
    exit
fi

. "$CONDA_PREFIX"/bin/activate # Source Conda base
conda activate teraseq

i="hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
sdir="$SAMPLE_DIR/$i"

mkdir -p "$sdir"/align || true # dliu
mkdir "$sdir"/fastq || true # dliu
mkdir "$sdir"/db || true # dliu

echo ">>> MERGE GENOME BAMS <<<"

inbam="reads.1.sanitize.toGenome.sorted.bam"

samtools merge -c -p -@ $threads "$sdir"/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/align/$inbam # dliu -c -p to prevent random suffixes
samtools index "$sdir"/align/$inbam &

echo ">>> MERGE TRANSCRIPTOME BAMS <<<"

inbam="reads.1.sanitize.noribo.toTranscriptome.sorted.bam"

samtools merge -c -p -@ $threads "$sdir"/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/align/$inbam # dliu -c -p to prevent random suffixes
samtools index "$sdir"/align/$inbam &
wait

echo ">>> MERGE W/O ADAPTER LISTS - REL5 <<<"

inlist="reads.1.sanitize.w_rel5.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > "$sdir"/fastq/$inlist

inlist="reads.1.sanitize.wo_rel5.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > "$sdir"/fastq/$inlist

echo ">>> MERGE W/O ADAPTER LISTS - REL3 <<<"

inlist="reads.1.sanitize.w_rel3.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > "$sdir"/fastq/$inlist

inlist="reads.1.sanitize.wo_rel3.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > "$sdir"/fastq/$inlist

echo ">>> MERGE DB FILES <<<"
# Since we have just a few reads we will merge them. They are more biological than technical replicates

samples="hsa.dRNASeq.HeLa.total.REL5.long.REL3.4 hsa.dRNASeq.HeLa.total.REL5.long.REL3.5 hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"

firstsample=""
for i in $samples; do
    firstsample="$firstsample$i"
    break
done

## Ttranscript table
for i in $samples; do
    echo " Working for" "$i"

    sqlite3 "$i"/db/sqlite.db ".dump transcr" > "$sdir"/db/"$i".sqlite.transcr.sql & #  make sqldump
done
wait

# Merge all sql dumps
# Use one of the samples to get the header
line_num=$(head -100 "$sdir"/db/"$firstsample".sqlite.transcr.sql | grep -nP "^CREATE TABLE transcr" | cut -d":" -f1) # Get line number/size of header
head -"${line_num}" "$sdir"/db/"$firstsample".sqlite.transcr.sql > "$sdir"/db/sqlite.transcr.sql # Get the header from ultimate 6

# Check the INDEX & COMMIT of the db and in case we are missing it append it manually
end_line=$(cat "$sdir"/db/"$firstsample".sqlite.transcr.sql | grep -nP "^CREATE INDEX transcr" | cut -d":" -f1) # Get line number/size of tail
if [ -z "$end_line" ]; then
    printf "CREATE INDEX transcr_loc ON transcr (rname, start);\nCOMMIT;\n" > "$sdir"/db/sqlite.transcr.sql.tail.tmp # We can just add it manually
    # dliu change echo -e to printf for POSIX compliance. trailing \n added.
else
    tail -n +"$end_line" "$sdir"/db/"$firstsample".sqlite.transcr.sql > "$sdir"/db/sqlite.transcr.sql.tail.tmp # Temporarily get the tail to append to the main db from ultimate 6
fi

add_num=200000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library
for i in $samples; do
    echo "$i"
#    add_num=$(( $(od -vAn -N2 -tu2 /dev/random) ))0000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library

    grep -P "^INSERT" "$sdir"/db/"$i".sqlite.transcr.sql | sed "s/INSERT INTO transcr VALUES(/INSERT INTO transcr VALUES($add_num/g" \
        >> "$sdir"/db/sqlite.transcr.sql # make unique id by adding ten-milions otherwise we get an error about not-unique id; keep it number makes it easier

    rm "$sdir"/db/"$i".sqlite.transcr.sql
    add_num=$((add_num + 200000000))
done

cat "$sdir"/db/sqlite.transcr.sql.tail.tmp >> "$sdir"/db/sqlite.transcr.sql && rm "$sdir"/db/sqlite.transcr.sql.tail.tmp # Append the tail and remove the temp

## genome table
for i in $samples; do
    echo " Working for" "$i"

    sqlite3 "$i"/db/sqlite.db ".dump genome" > "$sdir"/db/"$i".sqlite.genome.sql & #  make sqldump
done
wait

# Merge all sql dumps
# use one of the samples to get the header
line_num=$(head -100 "$sdir"/db/"$firstsample".sqlite.genome.sql | grep -nP "^CREATE TABLE genome" | cut -d":" -f1) # Get line number/size of header
head -"$line_num" "$sdir"/db/"$firstsample".sqlite.genome.sql > "$sdir"/db/sqlite.genome.sql # Get the header from ultimate 6

# Check the INDEX & COMMIT of the db and in case we are missing it append it manually
end_line=$(cat "$sdir"/db/"$firstsample".sqlite.genome.sql | grep -nP "^CREATE INDEX genome" | cut -d":" -f1) # Get line number/size of tail
if [ -z "$end_line" ]; then
    printf "CREATE INDEX genome_loc ON genome (rname, start);\nCOMMIT;\n" > "$sdir"/db/sqlite.genome.sql.tail.tmp # We can just add it manually
    # dliu change echo -e to printf for POSIX compliance. trailing \n added.
else
    tail -n +"$end_line" "$sdir"/db/"$firstsample".sqlite.genome.sql > "$sdir"/db/sqlite.genome.sql.tail.tmp # Temporarily get the tail to append to the main db from ultimate 6
fi

add_num=800000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library
for i in $samples; do
    echo "$i"
#    add_num=$(( $(od -vAn -N2 -tu2 /dev/random) ))0000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library

    # db=$i.sqlite.genome.sql # dliu unused

    grep -P "^INSERT" "$sdir"/db/"$i".sqlite.genome.sql | sed "s/INSERT INTO genome VALUES(/INSERT INTO genome VALUES($add_num/g" \
        >> "$sdir"/db/sqlite.genome.sql # make unique id by adding ten-milions otherwise we get an error about not-unique id; keep it number makes it easier
    rm "$sdir"/db/"$i".sqlite.genome.sql
    add_num=$((add_num + 200000000))
done

cat "$sdir"/db/sqlite.genome.sql.tail.tmp >> "$sdir"/db/sqlite.genome.sql && rm "$sdir"/db/sqlite.genome.sql.tail.tmp # Append the tail and remove the temp

# Make the "final" db
cat "$sdir"/db/sqlite.genome.sql "$sdir"/db/sqlite.transcr.sql | sqlite3 "$sdir"/db/sqlite.db && rm "$sdir"/db/sqlite.genome.sql "$sdir"/db/sqlite.transcr.sql

conda deactivate # dliu add deactivate to match unmatched activate

# Note: If you have the fast5 files you can continue with the analysis. Check https://github.com/mourelatos-lab/TERA-Seq_manuscript/samples/README.md for the location where to download them.

# echo ">>> MERGE NANOPOLISH POLYA <<<"
#
# inlist="reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab"
#
# if [ `ls -1 -a hsa.dRNASeq.HeLa.total.REL5.long.REL3.[4,5,6]/align/${inlist} 2>/dev/null | wc -l` -eq 3 ]; then
#     . $INSTALL/perl-virtualenv/teraseq/bin/activate
#
#     table-cat \
#         hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/align/$inlist \
#         hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/align/$inlist \
#         hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/align/$inlist \
#         > $sdir/align/$inlist
#
#     echo ">>> ANNOTATE WITH POLYA <<<"
#     ### Make list of polyA reads (>=0) to append them to db
#     # From Nanopolish:
#     #    qc_tag is an additional flag used to indicate the validity of the estimate. Generally speaking, you should only use rows of the output file with
#     #    this value set to PASS; all other rows with (e.g.) the qc_tag set to SUFFCLIP, ADAPTER, etc. display signs of irregularity that indicate that we believe the estimate to be unreliable.
#     inlist="reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab"
#
#     tail -n+2 $sdir/align/$inlist | egrep -w "PASS|NOREGION" \
#         | cut -f1,9 > $sdir/align/${inlist%.tab}.filt.tab
#
#     annotate-sqlite-with-file --db_col_add "polya" --db_col_bind "qname" \
#         --db_tables "genome" --database "$sdir/db/sqlite.db" --ifile "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" \
#         --round
#
#     annotate-sqlite-with-file --db_col_add "polya" --db_col_bind "qname" \
#         --db_tables "transcr" --database "$sdir/db/sqlite.db" --ifile "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" \
#         --round
# else
#     echo "It seems that $inlist doesn't exist for at least one of the samples. Please check you created hsa.dRNASeq.HeLa.total.REL5.long.REL3.[4,5,6]/fast5, downloaded and uncompressed fast5 tar.gz and placed the files in hsa.dRNASeq.HeLa.total.REL5.long.REL3.[4,5,6]/fast5 directories."
# fi

echo ">>> ALL DONE <<<"
