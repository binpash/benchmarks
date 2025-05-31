#!/usr/bin/env bash
set -x
set -e

TOP="$(git rev-parse --show-toplevel)"

echo ">>> SETUP DATA <<<"
#
# Prepare references and annotations
#
. "$TOP/bio/scripts/PARAMS.sh"

echo ">>> MAKE HG38 REFERENCES <<<"

assembly="hg38"

cd "$DATA_DIR" || exit 1

echo " >>> GET SILVA rRNA DATABASE <<<"

cd silva/ || exit 1

# Merge all ribosomal sequences together eliminating duplicates
zcat \
    SILVA_132_LSURef_tax_silva_trunc.fasta.gz \
    SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz \
    | fasta-rm-dup-seqs \
    | fasta-unique-names \
    | sed '/^[^>]/ y/uU/tT/' \
    > ribosomal.fa

# Create simple BED file where stop position corresponds to sequence length.
fasta-sanitize-header \
    --input ribosomal.fa \
    --delim : \
    --keep 0 \
    | fasta-to-sizes-bed \
        --name-suffix ribo \
        > ribosomal.bed

# Isolate human ribosomal genes and eliminate superfluous info from header.
grep -A 1 --no-group-separator "Homo sapiens" ribosomal.fa \
    | fasta-sanitize-header \
        --input - \
        --delim : \
        --keep 0 \
        > ribosomal.hsa.fa

echo ">>> GET GENOME AND REFERENCES <<<"
cd "$DATA_DIR/$assembly/" || exit 1

# Get chromosome sizes
samtools faidx genome/genome.fa
cut -f1-2 genome/genome.fa.fai > chrom.sizes

echo ">>> CLEAN ANNOTATION <<<"

# To keep protein-coding and polyA transcripts and delete transcripts that don't have either start or stop codons defined
cat ensembl_genes.gtf |
    clean-gtf-lines-polya --gtf - > genes-polya.gtf
ln -sf genes-polya.gtf genes.gtf

# Remove pseudogenes from annotation for toral samples
cat ensembl_genes.gtf |
    clean-gtf-lines-total --gtf - > genes-total.gtf

# Extract the genic elements (utr5, cds, utr3) from GFF and write them as BED.
gff-to-genic-elements-bed --input genes.gtf > genic_elements.bed || true

# Create separate files for each genic element.
grep -P ":utr5\t" genic_elements.bed > genic_elements.utr5.bed
grep -P ":utr3\t" genic_elements.bed > genic_elements.utr3.bed
grep -P ":cds\t" genic_elements.bed > genic_elements.cds.bed
grep -P ":ncrna\t" genic_elements.bed > genic_elements.ncrna.bed
grep -P ":mrna\t" genic_elements.bed > genic_elements.mrna.bed

# Create separate files for each genic element - total RNA - just copy the same thing as for polyA
ln -fs genic_elements.bed genic_elements-total.bed
ln -fs genic_elements.utr5.bed genic_elements-total.utr5.bed
ln -fs genic_elements.utr3.bed genic_elements-total.utr3.bed
ln -fs genic_elements.cds.bed genic_elements-total.cds.bed
ln -fs genic_elements.ncrna.bed genic_elements-total.ncrna.bed
ln -fs genic_elements.mrna.bed genic_elements-total.mrna.bed

echo ">>> MAKE PART OF MM10 REFERENCES <<<"
# We have to put this here because loading Conda messes up Perl libs
assembly="mm10"

cd "$DATA_DIR/silva/" || exit 1

grep -A 1 --no-group-separator "Mus musculus" ribosomal.fa \
    | fasta-sanitize-header \
        --input - \
        --delim : \
        --keep 0 \
        > ribosomal.mmu.fa

mkdir -p "$DATA_DIR/$assembly/genome"
cd "$DATA_DIR/$assembly/" || exit 1

echo ">>> CLEAN ANNOTATION <<<"

# To keep protein-coding and polyA transcripts and delete transcripts that don't have either start or stop codons defined
cat ensembl_genes.gtf \
    | clean-gtf-lines-polya --gtf - \
    > genes-polya.gtf
ln -sf genes-polya.gtf genes.gtf

echo ">>> RETURN TO HG38 REFERENCES <<<"
assembly="hg38"
cd "$DATA_DIR/$assembly" || exit 1

# Extract transcript sequences
gffread -w transcripts.fa -g genome/genome.fa genes.gtf # Poly(A)
gffread -w transcripts-total.fa -g genome/genome.fa genes-total.gtf # Total

echo ">>> ADD RRNA ANNOTATION <<<"
## Remove rRNA annotation from Ensembl and replace it with SILVA rRNA db annotation
# Make SILVA rRNA db to genome mapping

mkdir -p gmap-2019-09-12
gmap_build -d genome -D gmap-2019-09-12 genome/genome.fa # Make index
gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.hsa.fa -t $threads \
    --suboptimal-score=0.0 -f gff3_match_cdna |
    sed "s/\tcDNA_match\t/\texon\t/g" |
    sed "s/\tgenome\t/\tsilva\t/g" \
    > ribosomal.hsa.gmap.gff3 # Map ribosomal RNA to genome and get gff3 output to be added to the gene annotation
gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.hsa.fa -t $threads \
    --suboptimal-score=0.0 -S > ribosomal.hsa.gmap-summary.txt # Map ribosomal RNA to genome and get txt output for manual check
gff2gtf-gmap ribosomal.hsa.gmap.gff3 |
    sed "s/\tgmapidx\t/\tsilva\t/g" |
    sed "s/;$/; gene_biotype \"rRNA\"; transcript_biotype \"rRNA\";/g" |
    sort -k1,1 -k4,4n > ribosomal.hsa.gmap.gtf # Convert GMAP gff3 to gtf

gtf2bed6 ribosomal.hsa.gmap.gtf | cut -f1-6 > rRNA.bed

cat ensembl_genes.gtf | grep -v " \"rRNA\";" > ensembl_genes.gtf.tmp # Remove annotated rRNA from Ensembl but keep ribosomal - SILVA rRNA db doesn't annotate those
cat ensembl_genes.gtf.tmp ribosomal.hsa.gmap.gtf > ensembl_genes.gtf.tmp2
(grep "^#" ensembl_genes.gtf.tmp2; grep -v "^#" ensembl_genes.gtf.tmp2 | sort -k1,1 -k4,4n) > ensembl_genes.gtf # Add SILVA rRNA to Ensembl
rm ensembl_genes.gtf.tmp ensembl_genes.gtf.tmp2
gzip -c ensembl_genes.gtf > ensembl_genes.gtf.gz
ln -fs ensembl_genes.gtf.gz Homo_sapiens.GRCh38.91.gtf.gz

gffread -w ensembl-transcripts-wRibo.fa -g genome/genome.fa ensembl_genes.gtf # All the transcripts with rRNA

echo ">>> SUBSET PROTEIN-CODING TRANSCRIPTS <<<"

cat ensembl_genes.gtf | grep "transcript_biotype \"protein_coding\"" > ensembl_transcripts_protein_coding.gtf

echo ">>> MAKE MINIMAP2 INDEX <<<"

mkdir -p minimap2.17

ln -fs ../transcripts.fa minimap2.17/
ln -fs ../transcripts-total.fa minimap2.17/
ln -fs ../ensembl-transcripts-wRibo.fa  minimap2.17/
ln -fs ../genome/genome.fa minimap2.17/

minimap2 -k 12 -d minimap2.17/genome.k12.mmi minimap2.17/genome.fa
minimap2 -k 12 -d minimap2.17/transcripts.k12.mmi minimap2.17/transcripts.fa
minimap2 -k 12 -d minimap2.17/transcripts-total.k12.mmi minimap2.17/transcripts-total.fa
minimap2 -k 12 -d minimap2.17/ensembl-transcripts-wRibo.k12.mmi minimap2.17/ensembl-transcripts-wRibo.fa

# Convert hg38 to GRCh38 chromosome coding
cat GtRNAdb/hg38-tRNAs.bed |
    sed 's/^chrM/MT/g' |
    sed 's/^chr//g' |
    sed 's/chr1_KI270713v1_random/KI270713.1/g' |
    sort --parallel=$threads -T GtRNAdb/ -k1,1 -k2,2 |
    cut -f1-6 > tRNA.bed
cat tRNA.bed rRNA.bed > rRNA_tRNA.bed

echo ">>> GET CAGE SIGNALS <<<"

liftOver fantom5/HeLa.rep1.hg19.ctss_chr.bed.gz hg19ToHg38.over.chain.gz fantom5/HeLa.rep1.hg38.ctss_chr.bed fantom5/HeLa.rep1.hg19tohg38unmap.ctss_chr.bed &
liftOver fantom5/HeLa.rep2.hg19.ctss_chr.bed.gz hg19ToHg38.over.chain.gz fantom5/HeLa.rep2.hg38.ctss_chr.bed fantom5/HeLa.rep2.hg19tohg38unmap.ctss_chr.bed &
liftOver fantom5/HeLa.rep3.hg19.ctss_chr.bed.gz hg19ToHg38.over.chain.gz fantom5/HeLa.rep3.hg38.ctss_chr.bed fantom5/HeLa.rep3.hg19tohg38unmap.ctss_chr.bed &
wait

# Sort
RND=$RANDOM
cat fantom5/HeLa.rep1.hg38.ctss_chr.bed | sort --parallel=$threads -T . -k1,1 -k2,2n > tmp.$RND; mv tmp.$RND fantom5/HeLa.rep1.hg38.ctss_chr.bed
cat fantom5/HeLa.rep2.hg38.ctss_chr.bed | sort --parallel=$threads -T . -k1,1 -k2,2n > tmp.$RND; mv tmp.$RND fantom5/HeLa.rep2.hg38.ctss_chr.bed
cat fantom5/HeLa.rep3.hg38.ctss_chr.bed  | sort --parallel=$threads -T . -k1,1 -k2,2n > tmp.$RND; mv tmp.$RND fantom5/HeLa.rep3.hg38.ctss_chr.bed

# Convert from UCSC to Ensembl chromosomes naming
cat fantom5/HeLa.rep1.hg38.ctss_chr.bed | sed 's/^chrM/MT/g' | sed 's/^chr//g' | sed 's/14_GL000009v2_random/GL000009.2/g' \
    | sed 's/1_KI270706v1_random/KI270706.1/g' | sed 's/Un_KI270742v1/KI270742.1/g' | sort --parallel=$threads -T fantom5/ -k1,1 -k2,2n > fantom5/HeLa.rep1.hg38.ctss.bed
cat fantom5/HeLa.rep2.hg38.ctss_chr.bed | sed 's/^chrM/MT/g' | sed 's/^chr//g' | sed 's/14_GL000009v2_random/GL000009.2/g' \
    | sed 's/1_KI270706v1_random/KI270706.1/g' | sed 's/Un_KI270742v1/KI270742.1/g' | sort --parallel=$threads -T fantom5/ -k1,1 -k2,2n > fantom5/HeLa.rep2.hg38.ctss.bed
cat fantom5/HeLa.rep3.hg38.ctss_chr.bed | sed 's/^chrM/MT/g' | sed 's/^chr//g' | sed 's/14_GL000009v2_random/GL000009.2/g' \
    | sed 's/1_KI270706v1_random/KI270706.1/g' | sed 's/Un_KI270742v1/KI270742.1/g' | sort --parallel=$threads -T fantom5/ -k1,1 -k2,2n > fantom5/HeLa.rep3.hg38.ctss.bed

echo ">>> GET NET-CAGE SIGNALS <<<"

# Convert to hg19 to hg38
for gz in NET-CAGE/*.bed.gz; do
    base=${gz%.bed.gz}

    fifo=$(mktemp -u)
    mkfifo "$fifo"
    ( gunzip -c "$gz" >"$fifo" ) &

  #liftOver \
   # <(gunzip -c "$gz") \ 
    liftOver "$fifo" \
    hg19ToHg38.over.chain.gz \
    "${base}.hg38.ctss_chr.raw.bed" \
    "${base}.hg38.unmap.ctss_chr.bed"

  sort -k1,1 -k2,2n --parallel=$threads "${base}.hg38.ctss_chr.raw.bed" |
    sed \
      -e 's/^chrM/MT/' \
      -e 's/^chr//' \
      -e 's/14_GL000009v2_random/GL000009.2/' \
      -e 's/1_KI270706v1_random/KI270706.1/' \
      -e 's/Un_KI270742v1/KI270742.1/' |
    sort -k1,1 -k2,2n --parallel=$threads \
    > "${base}.hg38.ctss_chr.bed"

  rm "${base}.hg38.ctss_chr.raw.bed"
done

#cat meth/encodeCcreHela.bed | cut -f 10 | sort | uniq -c
#  20023 CTCF-only,CTCF-bound
#  31295 dELS
#   4169 dELS,CTCF-bound
#   2222 DNase-H3K4me3
#   1305 DNase-H3K4me3,CTCF-bound
#  35142 DNase-only
# 788464 Low-DNase
#  23551 pELS
#   4492 pELS,CTCF-bound
#  13231 PLS
#   2641 PLS,CTCF-bound
# Use mark "type" as name, not unique but easier to process
cat meth/encodeCcreHela.bed | cut -f 10 | tr ',' '\t' > meth/tmp
paste meth/encodeCcreHela.bed meth/tmp | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, $12, $5, $6}' \
    | substitute-in-column.py --table UCSC2ensembl.txt > meth/encodeCcreHela.genome.bed && rm meth/tmp # Convert UCSC chr to Ensembl
#cat meth/encodeCcreHela.genome.bed | cut -f 4 | sort | uniq -c
#  20023 CTCF-only
#  35464 dELS
#   3527 DNase-H3K4me3
#  35142 DNase-only
# 788464 Low-DNase
#  28043 pELS
#  15872 PLS

echo ">>> GET SIRV E2 REFERENCES <<<"

cd "$DATA_DIR/spikein/sirv/" || exit 1

sed -i 's/SIRVome_isoforms/SIRV/' SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta

# Get only annotation and fasta E2 from SIRV Set 1
# IMPORTANT: transcript fields in SIRV annotation often DON'T have correct strand
# This is how we can fix it
cat SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf | grep -P "\ttranscript\t|\texon\t" > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.tmp
fix-transcript-field-gtf SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.tmp SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.fixed.gtf && rm SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.tmp
ln -sf SIRVome_isoforms_C_170612a.fixed.gtf SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf
cat SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf | awk -v var="exon" 'BEGIN {FS="\t";OFS="\t"} {if ($3==var) {print $1, $4-1,$5, ".", ".", $7, $3, $9}}' \
    | sort -k1,1 -k2,2n > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.bed
gtf2bed12 SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.bed12
ln -sf SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf SIRV_Set1_Sequences_170612a/SIRV_isoforms_multi-fasta-annotation_C_170612a.E2.gtf

gffread -w SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa \
    -g SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf # Poly(A)

samtools faidx SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa
cat SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa.fai | cut -f1-2 \
    | sed '1i transcript_id\tlength' > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.length.tab

mkdir -p SIRV_Set1_Sequences_170612a/minimap2.17
ln -sf ../SIRVome_isoforms_170612a.fasta SIRV_Set1_Sequences_170612a/minimap2.17/
ln -sf ../SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa SIRV_Set1_Sequences_170612a/minimap2.17/

minimap2 -k 12 -d SIRV_Set1_Sequences_170612a/minimap2.17/transcripts_sirv1.E2.k12.mmi SIRV_Set1_Sequences_170612a/minimap2.17/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa
minimap2 -k 12 -d SIRV_Set1_Sequences_170612a/minimap2.17/genome.k12.mmi SIRV_Set1_Sequences_170612a/minimap2.17/SIRVome_isoforms_170612a.fasta

echo ">>> MAKE STAR INDEX <<<"
# Build STAR index on the genome without gene annotation
# Note: For human, STAR will require ~33 GB RAM unless you change --genomeSAsparseD settings

cd "$DATA_DIR/$assembly/" || exit 1

mkdir -p STAR-2.7.2b
STAR \
    --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir STAR-2.7.2b/ \
    --genomeFastaFiles genome/genome.fa \
    --genomeSAsparseD 3 # add this if you need to save RAM requirements; you can also increase the value to 3

# Build STAR index on the genome with gene annotation
mkdir -p STAR-2.7.2b-annot
STAR \
    --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir STAR-2.7.2b-annot/ \
    --genomeFastaFiles genome/genome.fa \
    --sjdbGTFfile ensembl_genes.gtf \
    --sjdbOverhang 100 \
    --genomeSAsparseD 3 # add this if you need to save RAM requirements; you can also increase the value to 3

echo ">>> MAKE MM10 REFERENCES <<<"

assembly="mm10"

cd "$DATA_DIR/$assembly/" || exit 1

echo ">>> ADD RRNA ANNOTATION <<<"
## Remove rRNA annotation from Ensembl and replace it with SILVA rRNA db annotation
# Make SILVA rRNA db to genome mapping

mkdir -p gmap-2019-09-12
gmap_build -d genome -D gmap-2019-09-12 genome/genome.fa # Make index
gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.mmu.fa -t $threads \
    --suboptimal-score=0.0 -f gff3_match_cdna | sed "s/\tcDNA_match\t/\texon\t/g" | sed "s/\tgenome\t/\tsilva\t/g" \
    > ribosomal.mmu.gmap.gff3 # Map ribosomal RNA to genome and get gff3 output to be added to the gene annotation
gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.mmu.fa -t $threads \
    --suboptimal-score=0.0 -S > ribosomal.mmu.gmap-summary.txt # Map ribosomal RNA to genome and get txt output for manual check
gff2gtf-gmap ribosomal.mmu.gmap.gff3 | sed "s/\tgmapidx\t/\tsilva\t/g" | sed "s/;$/; gene_biotype \"rRNA\"; transcript_biotype \"rRNA\";/g" \
    | sort -k1,1 -k4,4n > ribosomal.mmu.gmap.gtf # Convert GMAP gff3 to gtf

cat ensembl_genes.gtf | grep -v " \"rRNA\";" > ensembl_genes.gtf.tmp # Remove annotated rRNA from Ensembl but keep ribosomal - SILVA rRNA db doesn't annotate those
cat ensembl_genes.gtf.tmp ribosomal.mmu.gmap.gtf > ensembl_genes.gtf.tmp2
(grep "^#" ensembl_genes.gtf.tmp2; grep -v "^#" ensembl_genes.gtf.tmp2 | sort -k1,1 -k4,4n) > ensembl_genes.gtf # Add SILVA rRNA to Ensembl
rm ensembl_genes.gtf.tmp ensembl_genes.gtf.tmp2
gzip -c ensembl_genes.gtf > ensembl_genes.gtf.gz

echo ">>> ADD SIRV TO THE REFERENCES <<<"

# Add SIRV references to the references and do additional indexes
cat genome/genome.fa "$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta" > genome/genome_sirv1.fa

cat genes.gtf "$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf" > genes_sirv1.gtf

gffread -w transcripts_sirv1.fa -g genome/genome_sirv1.fa genes_sirv1.gtf

samtools faidx transcripts_sirv1.fa
grep -w -i sirv transcripts_sirv1.fa.fai | cut -f1-2 | sed '1i transcript_id\tlength' > sirv1.length.tab

mkdir -p minimap2.17
ln -sf ../transcripts_sirv1.fa minimap2.17/
ln -sf ../genome/genome_sirv1.fa minimap2.17/
minimap2 -k 12 -d minimap2.17/transcripts_sirv1.k12.mmi minimap2.17/transcripts_sirv1.fa
minimap2 -k 12 -d minimap2.17/genome_sirv1.k12.mmi minimap2.17/genome_sirv1.fa

# Get transcripts + ribosomal + sirv transcripts
cat ensembl_genes.gtf "$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf" > ensembl_genes_sirv1.gtf
gffread -w ensembl-transcripts-wRibo_sirv1.fa -g genome/genome_sirv1.fa ensembl_genes_sirv1.gtf

ln -sf ../ensembl-transcripts-wRibo_sirv1.fa  minimap2.17/
minimap2 -k 12 -d minimap2.17/ensembl-transcripts-wRibo_sirv1.k12.mmi minimap2.17/ensembl-transcripts-wRibo_sirv1.fa

echo ">>> ALL DONE <<<"
