#!/usr/bin/env bash
# TODO: Move this as part of the script
#
# process_data.sh
#   - Merge and sanitize FASTAs/BEDs
#   - Clean annotations, extract elements
#   - Build GMAP, minimap2 & STAR indices
#   - Liftover & format CAGE/NET-CAGE
#   - Prepare spike-in GTF/FASTA
#
set -e
set -x

TOP="$(git rev-parse --show-toplevel)"
source "$TOP/teraseq/scripts/PARAMS.sh"   # sets DATA_DIR, threads, etc.

################################################################################
# 1) SILVA rRNA processing
cd "$DATA_DIR/silva"
[[ ! -f ribosomal.fa ]] && \
  zcat SILVA_132_*tax_silva_trunc.fasta.gz \
    | fasta-rm-dup-seqs \
    | fasta-unique-names \
    | sed '/^[^>]/ y/uU/tT/' \
    > ribosomal.fa

[[ ! -f ribosomal.bed ]] && \
  fasta-sanitize-header --input ribosomal.fa --delim : --keep 0 \
    | fasta-to-sizes-bed --name-suffix ribo \
    > ribosomal.bed

# extract human & mouse rRNA
grep -A1 --no-group-separator "Homo sapiens" ribosomal.fa \
  | fasta-sanitize-header --input - --delim : --keep 0 > ribosomal.hsa.fa
grep -A1 --no-group-separator "Mus musculus" ribosomal.fa \
  | fasta-sanitize-header --input - --delim : --keep 0 > ribosomal.mmu.fa

################################################################################
# 2) hg38 references & annotations
cd "$DATA_DIR/hg38"
# clean Ensembl GTF → polyA vs total
[[ ! -f genes-polya.gtf ]] && \
  clean-gtf-lines-polya --gtf ensembl_genes.gtf > genes-polya.gtf
ln -sf genes-polya.gtf genes.gtf

[[ ! -f genes-total.gtf ]] && \
  clean-gtf-lines-total --gtf ensembl_genes.gtf > genes-total.gtf

# genic elements
[[ ! -f genic_elements.bed ]] && \
  gff-to-genic-elements-bed --input genes.gtf > genic_elements.bed

for elt in utr5 utr3 cds ncrna mrna; do
  grep -P ":${elt}\t" genic_elements.bed > genic_elements.${elt}.bed
  ln -fs genic_elements.${elt}.bed genic_elements-total.${elt}.bed
done

# GMAP rRNA → genome
cd "$DATA_DIR/hg38"
if [[ ! -d gmap-2019-09-12 ]]; then
  mkdir -p gmap-2019-09-12
  gmap_build -d genome -D gmap-2019-09-12 genome/genome.fa
  gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.hsa.fa -t $threads \
      --suboptimal-score=0.0 -f gff3_match_cdna | \
    sed -e 's/\tcDNA_match\t/\texon\t/g' -e 's/\tgenome\t/\tsilva\t/g' \
    > ribosomal.hsa.gmap.gff3
  gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.hsa.fa -t $threads \
      --suboptimal-score=0.0 -S \
    > ribosomal.hsa.gmap-summary.txt
fi

gff2gtf-gmap ribosomal.hsa.gmap.gff3 \
  | sed -e 's/\tgmapidx\t/\tsilva\t/g' \
        -e 's/;$/; gene_biotype "rRNA"; transcript_biotype "rRNA";/g' \
  | sort -k1,1 -k4,4n \
  > ribosomal.hsa.gmap.gtf
gtf2bed6 ribosomal.hsa.gmap.gtf | cut -f1-6 > rRNA.bed

# merge Ensembl + SILVA rRNA
grep -v ' "rRNA";' ensembl_genes.gtf > tmp.gtf
cat tmp.gtf ribosomal.hsa.gmap.gtf \
  | (grep '^#' && grep -v '^#' | sort -k1,1 -k4,4n) \
  > ensembl_genes.with_rRNA.gtf
mv ensembl_genes.with_rRNA.gtf ensembl_genes.gtf

# transcripts w/ & w/o rRNA
gffread -w transcripts.fa -g genome/genome.fa genes.gtf
gffread -w transcripts-total.fa -g genome/genome.fa genes-total.gtf
gffread -w ensembl-transcripts-wRibo.fa -g genome/genome.fa ensembl_genes.gtf

################################################################################
# 3) minimap2 indices
cd "$DATA_DIR"
for idx in hg38 hg38 hg38 hg38 mm10 spikein; do :; done  # placeholder
cd "$DATA_DIR"
mkdir -p minimap2.17
for ref in \
    hg38/genome/genome.fa \
    hg38/transcripts.fa \
    hg38/transcripts-total.fa \
    hg38/ensembl-transcripts-wRibo.fa \
    mm10/genome/genome.fa \
    spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta \
    spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa; do
  base=$(basename "${ref%.*}")
  minimap2 -k12 -d minimap2.17/${base}.k12.mmi "$ref"
done

################################################################################
# 4) STAR indices
cd "$DATA_DIR/hg38"
mkdir -p STAR-2.7.2b STAR-2.7.2b-annot
STAR --runMode genomeGenerate --runThreadN $threads \
     --genomeDir STAR-2.7.2b/ \
     --genomeFastaFiles genome/genome.fa \
     --genomeSAsparseD 3

STAR --runMode genomeGenerate --runThreadN $threads \
     --genomeDir STAR-2.7.2b-annot/ \
     --genomeFastaFiles genome/genome.fa \
     --sjdbGTFfile ensembl_genes.gtf \
     --sjdbOverhang 100 \
     --genomeSAsparseD 3

################################################################################
# 5) Liftover & reformat CAGE / NET-CAGE
cd "$DATA_DIR"
for f in fantom5/*.hg19.ctss_chr.bed.gz; do
  base=$(basename "$f" .bed.gz)
  gunzip -c "$f" \
    | liftOver stdin hg19ToHg38.over.chain.gz "${base}.hg38.ctss_chr.bed" /dev/null
  sed -e 's/^chrM/MT/' -e 's/^chr//' \
      -e 's/14_GL000009v2_random/GL000009.2/' \
      -e 's/1_KI270706v1_random/KI270706.1/' \
      -e 's/Un_KI270742v1/KI270742.1/' \
    < "${base}.hg38.ctss_chr.bed" \
    | sort -k1,1 -k2,2n \
    > fantom5/${base}.hg38.ctss.bed
done

for gz in NET-CAGE/*.ctss.bed.gz; do
  base=${gz%.ctss.bed.gz}
  gunzip -c "$gz" \
    | liftOver stdin hg19ToHg38.over.chain.gz ${base}.hg38.ctss_chr.raw.bed /dev/null
  sed -e 's/^chrM/MT/' -e 's/^chr//' \
      -e 's/14_GL000009v2_random/GL000009.2/' \
      -e 's/1_KI270706v1_random/KI270706.1/' \
      -e 's/Un_KI270742v1/KI270742.1/' \
    < ${base}.hg38.ctss_chr.raw.bed \
    | sort -k1,1 -k2,2n \
    > ${base}.hg38.ctss.bed
  rm ${base}.hg38.ctss_chr.raw.bed
done

################################################################################
# 6) Process ENCODE cCREs
cd "$DATA_DIR/meth"
cut -f10 encodeCcreHela.bed | tr ',' '\t' > tmp.types
paste encodeCcreHela.bed tmp.types \
  | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$12,$5,$6}' \
  | substitute-in-column.py --table ../UCSC2ensembl.txt \
  > encodeCcreHela.genome.bed
rm tmp.types

################################################################################
# 7) SPIKE-IN GTF/FASTA processing
cd "$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a"
grep -P "\t(transcript|exon)\t" SIRVome_isoforms_C_170612a.gtf \
  > tmp.gtf
fix-transcript-field-gtf tmp.gtf SIRVome_isoforms_C_170612a.fixed.gtf
mv SIRVome_isoforms_C_170612a.fixed.gtf SIRVome_isoforms_C_170612a.E2.gtf

gtf2bed12 SIRVome_isoforms_C_170612a.E2.gtf \
  > SIRVome_isoforms_C_170612a.E2.bed12

gffread -w SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa \
  -g SIRVome_isoforms_170612a.fasta \
  SIRVome_isoforms_C_170612a.E2.gtf

samtools faidx SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa
cut -f1-2 SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa.fai \
  | sed '1i transcript_id\tlength' \
  > SIRVome_isoforms_170612a.transcripts_sirv1.E2.length.tab

echo ">>> PROCESSING COMPLETE <<<"
