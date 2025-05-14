#!/bin/sh
#
# Default params for most of the analyses
# Can be overriden in individual scripts/steps
#

TERASEQ_ROOT="/root/TERA-Seq_manuscript"

CUR_DIR=$(pwd) # Current directory; useful for analyses
#DIR="$( cd "$( dirname "$0" )" && pwd )"
DIR="$( cd "$TERASEQ_ROOT" && pwd )"
INSTALL="$DIR/tools"
SAMPLE_DIR="$DIR/samples"
DATA_DIR="$DIR/data"
RES_DIR="results"

# Define common database filters
MATE1="((flag & 1) == 0 OR (flag & 64) == 64)"
MATE2="((flag & 1) == 0 OR (flag & 128) == 128)"
PRIMARY="((flag & 256) == 0)"
FORWARD="((flag & 16 ) == 0)"
REVERSE="((flag & 16) == 16)"
UNAMBIG="(best_hit == 1 AND number_of_best_hits == 1)"
UNIQUE="(number_of_mappings == 1)"
TRANSCRIPT="(transcript IS NOT NULL)"
CODINGTRANSCRIPT="(coding_transcript IS NOT NULL AND noncoding_transcript IS NULL)"
EXON="(exon IS NOT NULL)"
UTR5="(utr5 IS NOT NULL)"
CDS="(cds IS NOT NULL)"
UTR3="(utr3 IS NOT NULL)"
FEWMAPS="(number_of_mappings < 10)"

# Add utils directory to the path
export PATH="$DIR/tools/utils:${PATH}"

# common default assembly and number of threads
assembly="hg38"
threads=6
CONDA_PREFIX="/root/miniconda3"
