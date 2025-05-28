#!/usr/bin/env bash
set -x
set -e

TOP="$(git rev-parse --show-toplevel)"

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=full ;; # small uses a subset of full inputs
    --min) size=min ;;
    esac
done

export SIZE="$size" # for PARAMS.sh

SAMPLE_DIR="$TOP/teraseq/inputs/$size"
URL='https://atlas.cs.brown.edu/data'
samples="hsa.dRNASeq.HeLa.polyA.1 hsa.dRNASeq.HeLa.polyA.REL5.1 hsa.dRNASeq.HeLa.polyA.PNK.REL5.1"

base_url="$URL/teraseq/$size"

echo ">>> MAKE DIRECTORY STRUCTURE <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"
done

echo ">>> CHECK FASTQ <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for $i"

    if [ -f "$sdir/fastq/reads.1.fastq.gz" ]; then
        echo "$sdir/fastq/reads.1.fastq.gz is present, continue."
    else
        echo "$sdir/fastq/reads.1.fastq.gz does not exist, trying to download."
        download="$base_url/$i/fastq/reads.1.fastq.gz"
        mkdir -p "$sdir/fastq"
        curl "$download" > "$sdir/fastq/reads.1.fastq.gz"
    fi
done

echo ">>> SETUP DATA <<<"
#
# Prepare references and annotations
#
source "$TOP/teraseq/scripts/PARAMS.sh"

threads=8

####################################################################################################
echo ">>> MAKE HG38 REFERENCES <<<"

assembly="hg38"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR" || exit 1

echo " >>> GET SILVA rRNA DATABASE <<<"
# Download ribosomal sequences
# Nr99 "version" clusters highly (99%) similar sequences. From the Silva - "Ref NR 99 (Web database & ARB file), a 99% identity criterion to remove highly identical sequences using the  UCLUST tool was applied."
# Difference between Parc and Ref part of the database - "SILVA Parc and SILVA Ref. The Parc datasets comprise the entire SILVA databases for the respective gene, whereas the Ref datasets represent a subset of the Parc comprising only high-quality nearly full-length sequences."

mkdir -p silva
cd silva/ || exit 1

[ ! -f SILVA_132_LSURef_tax_silva_trunc.fasta.gz ] && \
    wget $base_url/data/SILVA_132_LSURef_tax_silva_trunc.fasta.gz

[ ! -f SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz ] && \
    wget $base_url/data/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz

echo ">>> GET GENOME AND REFERENCES <<<"

# Download and prepare the genome FASTA file from Ensembl
mkdir -p "$DATA_DIR/$assembly/genome"
cd "$DATA_DIR/$assembly/" || exit 1

[ ! -f genome/genome.fa ] && \
    wget -O- $base_url/data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz |
    gunzip -c |
    clean-genome-headers --fasta - > genome/genome.fa

# Download Ensembl annotation
[ ! -f ensembl_genes.gtf ] && \
    wget -O- $base_url/data/Homo_sapiens.GRCh38.91.gtf.gz |
    gunzip -c > ensembl_genes.gtf
ln -sf ensembl_genes.gtf Homo_sapiens.GRCh38.91.gtf

echo ">>> MAKE PART OF MM10 REFERENCES <<<"
assembly="mm10"

cd "$DATA_DIR/silva/" || exit 1

mkdir -p "$DATA_DIR/$assembly/genome"
cd "$DATA_DIR/$assembly/" || exit 1

[ ! -f genome/genome.fa ] && \
    wget -O- $base_url/data/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz |
    gunzip -c |
    clean-genome-headers --fasta - > genome/genome.fa

[ ! -f ensembl_genes.gtf ] && \
    wget -O- $base_url/data/Mus_musculus.GRCm38.97.gtf.gz |
    gunzip -c > ensembl_genes.gtf
ln -fs ensembl_genes.gtf Mus_musculus.GRCm38.97.gtf

echo ">>> RETURN TO HG38 REFERENCES <<<"

assembly="hg38"
cd "$DATA_DIR/$assembly" || exit 1

echo ">>> MAKE TRNA AND RRNA BED <<<"
# We can use this for cleaning of the bam files from tRNA and rRNA
mkdir -p GtRNAdb
[ ! -f GtRNAdb/hg38-tRNAs.tar.gz ] && \
    wget "$base_url/data/hg38-tRNAs.tar.gz" -O GtRNAdb/hg38-tRNAs.tar.gz
tar -xvzf GtRNAdb/hg38-tRNAs.tar.gz -C GtRNAdb

echo ">>> GET POLYA DATABASE <<<"

# PolyASite v2.0 (released 2019-08-13)
# IMPORTANT: PolyASite v2.0 is in GRCh38-96 Ensembl coordinates
mkdir -p polyasite-2.0
[ ! -f polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz ] && \
    wget $base_url/data/atlas.clusters.hg38.2-0.bed.gz -O polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz
gunzip -f polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz

echo ">>> GET CAGE SIGNALS <<<"
# CAGE data directly from FANTOM5 and convert to hg38 (hg19 in the database)
# Get FANTOM5 HeLa only
if [ ! -d fantom5 ]; then
    mkdir fantom5
    wget $base_url/data/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep1.CNhs12325.10815-111B5.hg19.ctss.bed.gz -O fantom5/HeLa.rep1.hg19.ctss_chr.bed.gz
    wget $base_url/data/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep2.CNhs12326.10816-111B6.hg19.ctss.bed.gz -O fantom5/HeLa.rep2.hg19.ctss_chr.bed.gz
    wget $base_url/data/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep3.CNhs12327.10817-111B7.hg19.ctss.bed.gz -O fantom5/HeLa.rep3.hg19.ctss_chr.bed.gz
fi

## Download required files for liftover from hg19 to hg38
[ ! -f hg19ToHg38.over.chain.gz ] && \
    wget $base_url/data/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz

echo ">>> GET NET-CAGE SIGNALS <<<"
# https://www.nature.com/articles/s41588-019-0485-9; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118075
if [ ! -d NET-CAGE ]; then
    mkdir NET-CAGE
    wget $base_url/data/NET-CAGE/Rep1-HeLaS3-NETCAGE-0_5M_CAC.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-0_5M_CAC.ctss.bed.gz
    wget $base_url/data/NET-CAGE/Rep1-HeLaS3-NETCAGE-1M_AGT.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-1M_AGT.ctss.bed.gz
    wget $base_url/data/NET-CAGE/Rep1-HeLaS3-NETCAGE-2M_GCG.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-2M_GCG.ctss.bed.gz
    wget $base_url/data/NET-CAGE/Rep2-HeLaS3-NETCAGE-0_5M_TAC.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-0_5M_TAC.ctss.bed.gz
    wget $base_url/data/NET-CAGE/Rep2-HeLaS3-NETCAGE-1M_ACG.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-1M_ACG.ctss.bed.gz
    wget $base_url/data/NET-CAGE/Rep2-HeLaS3-NETCAGE-2M_GCT.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-2M_GCT.ctss.bed.gz
fi

# Get conversion of UCSC->Ensembl
[ ! -f UCSC2ensembl.txt ] && \
    wget $base_url/data/UCSC2ensembl.txt -O UCSC2ensembl.txt

# Get cis-regions from ENCODE SEARCH https://screen.wenglab.org/
if [ ! -d meth ]; then
    mkdir meth
    wget $base_url/data/meth/encodeCcreHela.bed -O meth/encodeCcreHela.bed
fi

echo ">>> GET SIRV E2 REFERENCES <<<"
# Download Lexogen information
mkdir -p "$DATA_DIR/spikein/sirv"
cd "$DATA_DIR/spikein/sirv/" || exit 1

wget $URL/teraseq/TERA-Seq_manuscript/data/SIRV_Set1_Sequences_170612a.tar
tar -xvf SIRV_Set1_Sequences_170612a.tar

echo ">>> MAKE STAR INDEX <<<"
cd "$DATA_DIR/$assembly/" || exit 1

echo ">>> MAKE MM10 REFERENCES <<<"

assembly="mm10"

cd "$DATA_DIR/$assembly/" || exit 1

echo ">>> ALL DONE <<<"
