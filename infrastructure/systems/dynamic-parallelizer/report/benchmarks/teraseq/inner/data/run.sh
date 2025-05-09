#!/bin/bash
#
# Prepare references and annotations
#

source ../PARAMS.sh

threads=8

####################################################################################################
echo ">>> MAKE HG38 REFERENCES <<<"

assembly="hg38"

echo " >>> GET SILVA rRNA DATABASE <<<"
# Download ribosomal sequences
# Nr99 "version" clusters highly (99%) similar sequences. From the Silva - "Ref NR 99 (Web database & ARB file), a 99% identity criterion to remove highly identical sequences using the  UCLUST tool was applied."
# Difference between Parc and Ref part of the database - "SILVA Parc and SILVA Ref. The Parc datasets comprise the entire SILVA databases for the respective gene, whereas the Ref datasets represent a subset of the Parc comprising only high-quality nearly full-length sequences."

mkdir silva
cd silva/

wget \
    https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz \
    https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz

# Merge all ribosomal sequences together eliminating duplicates
source $INSTALL/perl-virtualenv/teraseq/bin/activate

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
# Download and prepare the genome FASTA file from Ensembl
mkdir -p $DATA_DIR/$assembly/genome
cd $DATA_DIR/$assembly/

wget -qO- ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
    | gunzip -c \
    | clean-genome-headers --fasta - \
    > genome/genome.fa

# Get chromosome sizes
$CONDA_PREFIX/envs/teraseq/bin/samtools faidx genome/genome.fa
cut -f1-2 genome/genome.fa.fai > chrom.sizes

# Download Ensembl annotation
wget -qO- ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz \
    | gunzip -c \
    > ensembl_genes.gtf
ln -s ensembl_genes.gtf Homo_sapiens.GRCh38.91.gtf

echo ">>> CLEAN ANNOTATION <<<"

# To keep protein-coding and polyA transcripts and delete transcripts that don't have either start or stop codons defined
cat ensembl_genes.gtf \
    | clean-gtf-lines-polya --gtf - \
    > genes-polya.gtf
ln -s genes-polya.gtf genes.gtf

# Remove pseudogenes from annotation for toral samples
cat ensembl_genes.gtf \
    | clean-gtf-lines-total --gtf - \
    > genes-total.gtf

# Extract the genic elements (utr5, cds, utr3) from GFF and write them as BED.
gff-to-genic-elements-bed \
    --input genes.gtf \
    > genic_elements.bed

# Create separate files for each genic element.
grep -P ":utr5\t" genic_elements.bed > genic_elements.utr5.bed
grep -P ":utr3\t" genic_elements.bed > genic_elements.utr3.bed
grep -P ":cds\t" genic_elements.bed > genic_elements.cds.bed
grep -P ":ncrna\t" genic_elements.bed > genic_elements.ncrna.bed
grep -P ":mrna\t" genic_elements.bed > genic_elements.mrna.bed

# Create separate files for each genic element - total RNA - just copy the same thing as for polyA
ln -s genic_elements.bed genic_elements-total.bed
ln -s genic_elements.utr5.bed genic_elements-total.utr5.bed
ln -s genic_elements.utr3.bed genic_elements-total.utr3.bed
ln -s genic_elements.cds.bed genic_elements-total.cds.bed
ln -s genic_elements.ncrna.bed genic_elements-total.ncrna.bed
ln -s genic_elements.mrna.bed genic_elements-total.mrna.bed

echo ">>> MAKE PART OF MM10 REFERENCES <<<"
# We have to put this here because loading Conda messes up Perl libs
assembly="mm10"

cd $DATA_DIR/silva/

grep -A 1 --no-group-separator "Mus musculus" ribosomal.fa \
    | fasta-sanitize-header \
        --input - \
        --delim : \
        --keep 0 \
        > ribosomal.mmu.fa

mkdir -p $DATA_DIR/$assembly/genome
cd $DATA_DIR/$assembly/

wget -qO- ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
    | gunzip -c \
    | clean-genome-headers --fasta - \
    > genome/genome.fa

wget -qO- ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz | gunzip -c > ensembl_genes.gtf
ln -s ensembl_genes.gtf Mus_musculus.GRCm38.97.gtf

echo ">>> CLEAN ANNOTATION <<<"

# To keep protein-coding and polyA transcripts and delete transcripts that don't have either start or stop codons defined
cat ensembl_genes.gtf \
    | clean-gtf-lines-polya --gtf - \
    > genes-polya.gtf
ln -s genes-polya.gtf genes.gtf

deactivate

echo ">>> RETURN TO HG38 REFERENCES <<<"

assembly="hg38"
cd $DATA_DIR/$assembly/

# Extract transcript sequences
if [ -z ${CONDA_PREFIX} ]; then
    echo "Variable \$CONDA_PREFIX is not set. Please make sure you specified if in PARAMS.sh."
    exit
fi

source $CONDA_PREFIX/bin/activate # Source Conda base
conda activate teraseq

gffread -w transcripts.fa -g genome/genome.fa genes.gtf # Poly(A)
gffread -w transcripts-total.fa -g genome/genome.fa genes-total.gtf # Total

echo ">>> ADD RRNA ANNOTATION <<<"
## Remove rRNA annotation from Ensembl and replace it with SILVA rRNA db annotation
# Make SILVA rRNA db to genome mapping

mkdir gmap-2019-09-12
gmap_build -d genome -D gmap-2019-09-12 genome/genome.fa # Make index
gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.hsa.fa -t $threads \
    --suboptimal-score=0.0 -f gff3_match_cdna | sed "s/\tcDNA_match\t/\texon\t/g" | sed "s/\tgenome\t/\tsilva\t/g" \
    > ribosomal.hsa.gmap.gff3 # Map ribosomal RNA to genome and get gff3 output to be added to the gene annotation
gmap -d genome -D gmap-2019-09-12 ../silva/ribosomal.hsa.fa -t $threads \
    --suboptimal-score=0.0 -S > ribosomal.hsa.gmap-summary.txt # Map ribosomal RNA to genome and get txt output for manual check
gff2gtf-gmap ribosomal.hsa.gmap.gff3 | sed "s/\tgmapidx\t/\tsilva\t/g" | sed "s/;$/; gene_biotype \"rRNA\"; transcript_biotype \"rRNA\";/g" \
    | sort -k1,1 -k4,4n > ribosomal.hsa.gmap.gtf # Convert GMAP gff3 to gtf

gtf2bed6 ribosomal.hsa.gmap.gtf | cut -f1-6 > rRNA.bed

cat ensembl_genes.gtf | grep -v " \"rRNA\";" > ensembl_genes.gtf.tmp # Remove annotated rRNA from Ensembl but keep ribosomal - SILVA rRNA db doesn't annotate those
cat ensembl_genes.gtf.tmp ribosomal.hsa.gmap.gtf > ensembl_genes.gtf.tmp2
(grep "^#" ensembl_genes.gtf.tmp2; grep -v "^#" ensembl_genes.gtf.tmp2 | sort -k1,1 -k4,4n) > ensembl_genes.gtf # Add SILVA rRNA to Ensembl
rm ensembl_genes.gtf.tmp ensembl_genes.gtf.tmp2
gzip -c ensembl_genes.gtf > ensembl_genes.gtf.gz
ln -s ensembl_genes.gtf.gz Homo_sapiens.GRCh38.91.gtf.gz

gffread -w ensembl-transcripts-wRibo.fa -g genome/genome.fa ensembl_genes.gtf # All the transcripts with rRNA

echo ">>> SUBSET PROTEIN-CODING TRANSCRIPTS <<<"

cat ensembl_genes.gtf | grep "transcript_biotype \"protein_coding\"" > ensembl_transcripts_protein_coding.gtf

echo ">>> MAKE MINIMAP2 INDEX <<<"

mkdir -p minimap2.17

ln -s ../transcripts.fa minimap2.17/
ln -s ../transcripts-total.fa minimap2.17/
ln -s ../ensembl-transcripts-wRibo.fa  minimap2.17/
ln -s ../genome/genome.fa minimap2.17/

minimap2 -k 12 -d minimap2.17/genome.k12.mmi minimap2.17/genome.fa &
minimap2 -k 12 -d minimap2.17/transcripts.k12.mmi minimap2.17/transcripts.fa &
minimap2 -k 12 -d minimap2.17/transcripts-total.k12.mmi minimap2.17/transcripts-total.fa &
minimap2 -k 12 -d minimap2.17/ensembl-transcripts-wRibo.k12.mmi minimap2.17/ensembl-transcripts-wRibo.fa &
wait

echo ">>> MAKE TRNA AND RRNA BED <<<"
# We can use this for cleaning of the bam files from tRNA and rRNA
mkdir GtRNAdb
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz -O GtRNAdb/hg38-tRNAs.tar.gz
tar -xvzf GtRNAdb/hg38-tRNAs.tar.gz -C GtRNAdb
# Convert hg38 to GRCh38 chromosome coding
cat GtRNAdb/hg38-tRNAs.bed | sed 's/^chrM/MT/g' | sed 's/^chr//g' | sed 's/chr1_KI270713v1_random/KI270713.1/g' \
    | sort --parallel=$threads -T GtRNAdb/ -k1,1 -k2,2 | cut -f1-6 > tRNA.bed

cat tRNA.bed rRNA.bed > rRNA_tRNA.bed

echo ">>> GET POLYA DATABASE <<<"

# PolyASite v2.0 (released 2019-08-13)
# IMPORTANT: PolyASite v2.0 is in GRCh38-96 Ensembl coordinates
mkdir -p polyasite-2.0
wget https://polyasite.unibas.ch/download/clusters/GRCh38-96/2-0/atlas.clusters.hg38.2-0.bed.gz -O polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz
gunzip polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz

echo ">>> GET CAGE SIGNALS <<<"
# CAGE data directly from FANTOM5 and convert to hg38 (hg19 in the database)
# Get FANTOM5 HeLa only
mkdir fantom5
#wget http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep1.CNhs12325.10815-111B5.hg19.ctss.bed.gz -O fantom5/HeLa.rep1.hg19.ctss_chr.bed.gz
#wget http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep2.CNhs12326.10816-111B6.hg19.ctss.bed.gz -O fantom5/HeLa.rep2.hg19.ctss_chr.bed.gz
#wget http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep3.CNhs12327.10817-111B7.hg19.ctss.bed.gz -O fantom5/HeLa.rep3.hg19.ctss_chr.bed.gz
wget https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep1.CNhs12325.10815-111B5.hg19.ctss.bed.gz -O fantom5/HeLa.rep1.hg19.ctss_chr.bed.gz
wget https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep2.CNhs12326.10816-111B6.hg19.ctss.bed.gz -O fantom5/HeLa.rep2.hg19.ctss_chr.bed.gz
wget https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep3.CNhs12327.10817-111B7.hg19.ctss.bed.gz -O fantom5/HeLa.rep3.hg19.ctss_chr.bed.gz

## Download required files for liftover from hg19 to hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz
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
# https://www.nature.com/articles/s41588-019-0485-9; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118075
mkdir NET-CAGE
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318225/suppl/GSM3318225_CNhi10918_biologicalRep1-HeLaS3-NETCAGE-0_5M_CAC.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-0_5M_CAC.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318226/suppl/GSM3318226_CNhi10918_biologicalRep1-HeLaS3-NETCAGE-1M_AGT.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-1M_AGT.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318227/suppl/GSM3318227_CNhi10918_biologicalRep1-HeLaS3-NETCAGE-2M_GCG.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-2M_GCG.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318228/suppl/GSM3318228_CNhi10918_biologicalRep2-HeLaS3-NETCAGE-0_5M_TAC.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-0_5M_TAC.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318229/suppl/GSM3318229_CNhi10918_biologicalRep2-HeLaS3-NETCAGE-1M_ACG.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-1M_ACG.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318230/suppl/GSM3318230_CNhi10918_biologicalRep2-HeLaS3-NETCAGE-2M_GCT.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-2M_GCT.ctss.bed.gz

# Convert to hg19 to hg38
for i in NET-CAGE/*.bed.gz; do
    liftOver $i hg19ToHg38.over.chain.gz ${i%.ctss*}.hg38.ctss_chr.bed ${i%.ctss*}.hg38.unmap.ctss_chr.bed

    cat ${i%.ctss*}.hg38.ctss_chr.bed | sort --parallel=$threads -T . -k1,1 -k2,2n > tmp.$RND
    # Convert from UCSC to Ensembl chromosomes naming
    cat tmp.$RND | sed 's/^chrM/MT/g' | sed 's/^chr//g' | sed 's/14_GL000009v2_random/GL000009.2/g' \
        | sed 's/1_KI270706v1_random/KI270706.1/g' | sed 's/Un_KI270742v1/KI270742.1/g' | sort --parallel=$threads -T . -k1,1 -k2,2n > ${i%.ctss*}.hg38.ctss_chr.bed
    rm tmp.$RND
done

# Get conversion of UCSC->Ensembl
wget https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_UCSC2ensembl.txt -O UCSC2ensembl.txt

# dliu api.wenglab.org not reachable. skipping
# # Get cis-regions from ENCODE SEARCH https://screen.wenglab.org/
# mkdir meth
# wget https://api.wenglab.org/screen_v13/fdownloads/Seven-Group/ENCFF977IGB_ENCFF489CIY_ENCFF194XTD_ENCFF836JPY.7group.bed -O meth/encodeCcreHela.bed
# #cat meth/encodeCcreHela.bed | cut -f 10 | sort | uniq -c
# #  20023 CTCF-only,CTCF-bound
# #  31295 dELS
# #   4169 dELS,CTCF-bound
# #   2222 DNase-H3K4me3
# #   1305 DNase-H3K4me3,CTCF-bound
# #  35142 DNase-only
# # 788464 Low-DNase
# #  23551 pELS
# #   4492 pELS,CTCF-bound
# #  13231 PLS
# #   2641 PLS,CTCF-bound
# # Use mark "type" as name, not unique but easier to process
# cat meth/encodeCcreHela.bed | cut -f 10 | tr ',' '\t' > meth/tmp
# paste meth/encodeCcreHela.bed meth/tmp | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, $12, $5, $6}' \
#     | substitute-in-column.py --table UCSC2ensembl.txt > meth/encodeCcreHela.genome.bed && rm meth/tmp # Convert UCSC chr to Ensembl
# #cat meth/encodeCcreHela.genome.bed | cut -f 4 | sort | uniq -c
# #  20023 CTCF-only
# #  35464 dELS
# #   3527 DNase-H3K4me3
# #  35142 DNase-only
# # 788464 Low-DNase
# #  28043 pELS
# #  15872 PLS

echo ">>> GET SIRV E2 REFERENCES <<<"
# Download Lexogen information
mkdir -p $DATA_DIR/spikein/sirv
cd $DATA_DIR/spikein/sirv/

#wget https://www.lexogen.com/wp-content/uploads/2018/08/SIRV_Set1_Sequences_170612a-ZIP.zip # The original download link - you would need to uncomment the following section and install `unzip` and `dos2unix`
# unzip $DATA_DIR/SIRV_Set1_Sequences_170612a-ZIP.zip -d .
# mv SIRV_Set1_Sequences_170612a\ \(ZIP\) SIRV_Set1_Sequences_170612a
# for i in SIRV_Set1_Sequences_170612a/*; do
#     dos2unix $i
# done
tar -xvf $DATA_DIR/SIRV_Set1_Sequences_170612a.tar # dliu prefix with DATA_DIR because file otherwise wouldn't be found otherwise
sed -i 's/SIRVome_isoforms/SIRV/' SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta

# Get only annotation and fasta E2 from SIRV Set 1
# IMPORTANT: transcript fields in SIRV annotation often DON'T have correct strand
# This is how we can fix it
cat SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf | grep -P "\ttranscript\t|\texon\t" > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.tmp
fix-transcript-field-gtf SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.tmp SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.fixed.gtf && rm SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.tmp
ln -s SIRVome_isoforms_C_170612a.fixed.gtf SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf
cat SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf | awk -v var="exon" 'BEGIN {FS="\t";OFS="\t"} {if ($3==var) {print $1, $4-1,$5, ".", ".", $7, $3, $9}}' \
    | sort -k1,1 -k2,2n > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.bed
gtf2bed12 SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.bed12
ln -s SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf SIRV_Set1_Sequences_170612a/SIRV_isoforms_multi-fasta-annotation_C_170612a.E2.gtf

gffread -w SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa \
    -g SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf # Poly(A)

samtools faidx SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa
cat SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa.fai | cut -f1-2 \
    | sed '1i transcript_id\tlength' > SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.length.tab

mkdir SIRV_Set1_Sequences_170612a/minimap2.17
ln -s ../SIRVome_isoforms_170612a.fasta SIRV_Set1_Sequences_170612a/minimap2.17/
ln -s ../SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa SIRV_Set1_Sequences_170612a/minimap2.17/

minimap2 -k 12 -d SIRV_Set1_Sequences_170612a/minimap2.17/transcripts_sirv1.E2.k12.mmi SIRV_Set1_Sequences_170612a/minimap2.17/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa &
minimap2 -k 12 -d SIRV_Set1_Sequences_170612a/minimap2.17/genome.k12.mmi SIRV_Set1_Sequences_170612a/minimap2.17/SIRVome_isoforms_170612a.fasta &
wait

echo ">>> MAKE STAR INDEX <<<"
# Build STAR index on the genome without gene annotation
# Note: For human, STAR will require ~33 GB RAM unless you change --genomeSAsparseD settings

cd $DATA_DIR/$assembly/

mkdir STAR-2.7.2b
STAR \
    --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir STAR-2.7.2b/ \
    --genomeFastaFiles genome/genome.fa
#    --genomeSAsparseD 2 # add this if you need to save RAM requirements; you can also increase the value to 3

# Build STAR index on the genome with gene annotation
mkdir STAR-2.7.2b-annot
STAR \
    --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir STAR-2.7.2b-annot/ \
    --genomeFastaFiles genome/genome.fa \
    --sjdbGTFfile ensembl_genes.gtf \
    --sjdbOverhang 100
#    --genomeSAsparseD 2 # add this if you need to save RAM requirements; you can also increase the value to 3

echo ">>> MAKE MM10 REFERENCES <<<"

assembly="mm10"

cd $DATA_DIR/$assembly/

echo ">>> ADD RRNA ANNOTATION <<<"
## Remove rRNA annotation from Ensembl and replace it with SILVA rRNA db annotation
# Make SILVA rRNA db to genome mapping

mkdir gmap-2019-09-12
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
cat genome/genome.fa $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.fasta > genome/genome_sirv1.fa

cat genes.gtf $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf > genes_sirv1.gtf

gffread -w transcripts_sirv1.fa -g genome/genome_sirv1.fa genes_sirv1.gtf

samtools faidx transcripts_sirv1.fa
grep -w -i sirv transcripts_sirv1.fa.fai | cut -f1-2 | sed '1i transcript_id\tlength' > sirv1.length.tab

mkdir minimap2.17
ln -s ../transcripts_sirv1.fa minimap2.17/
ln -s ../genome/genome_sirv1.fa minimap2.17/
minimap2 -k 12 -d minimap2.17/transcripts_sirv1.k12.mmi minimap2.17/transcripts_sirv1.fa &
minimap2 -k 12 -d minimap2.17/genome_sirv1.k12.mmi minimap2.17/genome_sirv1.fa &
wait

# Get transcripts + ribosomal + sirv transcripts
cat ensembl_genes.gtf $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf > ensembl_genes_sirv1.gtf
gffread -w ensembl-transcripts-wRibo_sirv1.fa -g genome/genome_sirv1.fa ensembl_genes_sirv1.gtf

ln -s ../ensembl-transcripts-wRibo_sirv1.fa  minimap2.17/
minimap2 -k 12 -d minimap2.17/ensembl-transcripts-wRibo_sirv1.k12.mmi minimap2.17/ensembl-transcripts-wRibo_sirv1.fa &
wait

echo ">>> ALL DONE <<<"
