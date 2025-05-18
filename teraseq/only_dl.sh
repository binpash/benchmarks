#!/usr/bin/env bash
set -x

TOP="$(git rev-parse --show-toplevel)"

SAMPLE_DIR="$TOP/teraseq/inputs"
URL='https://atlas.cs.brown.edu/data'
outdir="$TOP/teraseq/outputs"
benchmark_dir="$TOP/teraseq"

samples="hsa.dRNASeq.HeLa.polyA.1 hsa.dRNASeq.HeLa.polyA.REL5.1 hsa.dRNASeq.HeLa.polyA.PNK.REL5.1"

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=small ;;
    --min) size=min ;;
    esac
done

echo ">>> MAKE DIRECTORY STRUCTURE <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for" "$i"

    mkdir -p "$sdir/logfiles"
    mkdir "$sdir/align"
    mkdir "$sdir/db"
done

echo ">>> CHECK FASTQ <<<"

for i in $samples; do
    sdir=$SAMPLE_DIR/$i
    echo " Working for $i"

    if [ -f "$sdir/fastq/reads.1.fastq.gz" ]; then
        echo "$sdir/fastq/reads.1.fastq.gz is present, continue."
    else
        echo "$sdir/fastq/reads.1.fastq.gz does not exist, trying to download."
        download="$(grep download "$benchmark_dir/README.md" | grep "$i" | cut -d '|' -f 6 | cut -d '(' -f2  | sed 's/)//' | sed 's#https://##' | tr -d '[:space:]')"
        mkdir -p "$sdir/fastq"
        curl "$download" > "$sdir/fastq/reads.1.fastq.gz"
    fi
done

wget $URL/teraseq.new/TERA-Seq_manuscript/data/SIRV_Set1_Sequences_170612a.tar

wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz

wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 

wget ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz

wget ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/mus_musculus.grcm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz 

wget 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz' -O GtRNAdb/hg38-tRNAs.tar.gz

wget https://polyasite.unibas.ch/download/clusters/GRCh38-96/2-0/atlas.clusters.hg38.2-0.bed.gz -O polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz

wget https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep1.CNhs12325.10815-111B5.hg19.ctss.bed.gz -O fantom5/HeLa.rep1.hg19.ctss_chr.bed.gz
wget https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep2.CNhs12326.10816-111B6.hg19.ctss.bed.gz -O fantom5/HeLa.rep2.hg19.ctss_chr.bed.gz
wget https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/epitheloid%2520carcinoma%2520cell%2520line%253a%2520HelaS3%2520ENCODE%252c%2520biol_rep3.CNhs12327.10817-111B7.hg19.ctss.bed.gz -O fantom5/HeLa.rep3.hg19.ctss_chr.bed.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318225/suppl/GSM3318225_CNhi10918_biologicalRep1-HeLaS3-NETCAGE-0_5M_CAC.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-0_5M_CAC.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318226/suppl/GSM3318226_CNhi10918_biologicalRep1-HeLaS3-NETCAGE-1M_AGT.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-1M_AGT.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318227/suppl/GSM3318227_CNhi10918_biologicalRep1-HeLaS3-NETCAGE-2M_GCG.ctss.bed.gz -O NET-CAGE/Rep1-HeLaS3-NETCAGE-2M_GCG.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318228/suppl/GSM3318228_CNhi10918_biologicalRep2-HeLaS3-NETCAGE-0_5M_TAC.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-0_5M_TAC.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318229/suppl/GSM3318229_CNhi10918_biologicalRep2-HeLaS3-NETCAGE-1M_ACG.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-1M_ACG.ctss.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3318nnn/GSM3318230/suppl/GSM3318230_CNhi10918_biologicalRep2-HeLaS3-NETCAGE-2M_GCT.ctss.bed.gz -O NET-CAGE/Rep2-HeLaS3-NETCAGE-2M_GCT.ctss.bed.gz

wget https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_UCSC2ensembl.txt -O UCSC2ensembl.txt

wget https://api.wenglab.org/screen_v13/fdownloads/Seven-Group/ENCFF977IGB_ENCFF489CIY_ENCFF194XTD_ENCFF836JPY.7group.bed -O meth/encodeCcreHela.bed
