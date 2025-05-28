#!/usr/bin/env bash
set -e
set -x
#
# download_data.sh
#   - Create directories
#   - Download FASTQ samples
#   - Download SILVA databases, genomes, GTFs, spike-ins, CAGE/NET-CAGE, PolyASite, etc.
#
set -euo pipefail

TOP="$(git rev-parse --show-toplevel)"
size=full
for arg in "$@"; do
  case "$arg" in
    --small) size=small ;;
    --min)   size=min ;;
  esac
done

export SIZE="$size"            # used by PARAMS.sh downstream
SAMPLE_DIR="$TOP/teraseq/inputs/$SIZE"
URL='https://atlas.cs.brown.edu/data'
base_url="$URL/teraseq/$SIZE"
samples=(
  hsa.dRNASeq.HeLa.polyA.1
  hsa.dRNASeq.HeLa.polyA.REL5.1
  hsa.dRNASeq.HeLa.polyA.PNK.REL5.1
)

echo ">>> MAKING DIRECTORY STRUCTURE <<<"
mkdir -p "$SAMPLE_DIR"
for s in "${samples[@]}"; do
  echo " Creating sample dir: $s"
  mkdir -p "$SAMPLE_DIR/$s/fastq"
done

echo ">>> DOWNLOADING FASTQ <<<"
for s in "${samples[@]}"; do
  fq="$SAMPLE_DIR/$s/fastq/reads.1.fastq.gz"
  if [[ -f "$fq" ]]; then
    echo "  $fq already exists"
  else
    echo "  Downloading $s FASTQ..."
    curl -fSL "$base_url/$s/fastq/reads.1.fastq.gz" -o "$fq"
  fi
done

echo ">>> PREPARING DATA_DIR & PARAMS <<<"
source "$TOP/teraseq/scripts/PARAMS.sh"  # sets DATA_DIR, etc.
mkdir -p "$DATA_DIR"

echo ">>> DOWNLOADING SILVA rRNA DATABASE <<<"
mkdir -p "$DATA_DIR/silva" && cd "$DATA_DIR/silva"
for db in SILVA_132_LSURef_tax_silva_trunc.fasta.gz SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz; do
  if [[ ! -s "$db" ]]; then
    wget "$base_url/data/$db"
  else
    echo "  $db exists"
  fi
done

echo ">>> DOWNLOADING GENOME REFERENCES <<<"
for asm in hg38/mm10; do
  target="$DATA_DIR/$asm"
  mkdir -p "$target/genome"
  cd "$target"
  if [[ ! -s genome/genome.fa ]]; then
    wget -O- "$base_url/data/$( [ "$asm" == "hg38" ] \
        && echo Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
        || echo Mus_musculus.GRCm38.dna.primary_assembly.fa.gz)" \
      | gunzip -c \
      | clean-genome-headers --fasta - \
      > genome/genome.fa
  fi

  if [[ ! -s ensembl_genes.gtf ]]; then
    wget -O- "$base_url/data/$( [ "$asm" == "hg38" ] \
        && echo Homo_sapiens.GRCh38.91.gtf.gz \
        || echo Mus_musculus.GRCm38.97.gtf.gz)" \
      | gunzip -c \
      > ensembl_genes.gtf
  fi
done

echo ">>> DOWNLOADING PolyASite, CAGE, NET-CAGE, UCSC→Ensembl mappings <<<"
cd "$DATA_DIR" || exit 1
# PolyASite v2.0
mkdir -p "$DATA_DIR/polyasite-2.0"
wget -O polyasite-2.0/atlas.clusters.hg38.2-0.bed.gz \
     "$base_url/data/atlas.clusters.hg38.2-0.bed.gz"

# FANTOM5 HeLa CAGE
mkdir -p "$DATA_DIR/fantom5"
for rep in 1 2 3; do
  file="HeLa.rep${rep}.hg19.ctss_chr.bed.gz"
  if [[ ! -s fantom5/$file ]]; then
    wget -O fantom5/$file \
      "$base_url/data/epitheloid-carcinoma-cell-line--HelaS3-ENCODE-biol_rep${rep}.CNhs1232${4+rep}.1081${5+rep}-111B${4+rep}.hg19.ctss.bed.gz"
  fi
done
wget -O hg19ToHg38.over.chain.gz "$base_url/data/hg19ToHg38.over.chain.gz"

# NET-CAGE
mkdir -p "$DATA_DIR/NET-CAGE"
for f in \
  Rep1-HeLaS3-NETCAGE-0_5M_CAC \
  Rep1-HeLaS3-NETCAGE-1M_AGT \
  Rep1-HeLaS3-NETCAGE-2M_GCG \
  Rep2-HeLaS3-NETCAGE-0_5M_TAC \
  Rep2-HeLaS3-NETCAGE-1M_ACG \
  Rep2-HeLaS3-NETCAGE-2M_GCT; do
  gz="$f.ctss.bed.gz"
  if [[ ! -s NET-CAGE/$gz ]]; then
    wget -O NET-CAGE/$gz "$base_url/data/NET-CAGE/$gz"
  fi
done

# UCSC2ensembl mapping & ENCODE cCREs
wget -O "$DATA_DIR/UCSC2ensembl.txt" "$base_url/data/UCSC2ensembl.txt"
mkdir -p "$DATA_DIR/meth"
wget -O "$DATA_DIR/meth/encodeCcreHela.bed" "$base_url/data/meth/encodeCcreHela.bed"

echo ">>> DOWNLOADING SIRV E2 SPIKE-INS <<<"
mkdir -p "$DATA_DIR/spikein/sirv" && cd "$DATA_DIR/spikein/sirv"
if [[ ! -f SIRV_Set1_Sequences_170612a.tar ]]; then
  wget "$URL/teraseq/TERA-Seq_manuscript/data/SIRV_Set1_Sequences_170612a.tar"
  tar -xf SIRV_Set1_Sequences_170612a.tar
fi

echo ">>> ALL DOWNLOADS COMPLETE <<<"
