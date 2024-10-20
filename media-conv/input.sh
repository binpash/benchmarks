#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
input_dir="${eval_dir}/input"

DATA_LINK=https://atlas-group.cs.brown.edu/data/wav.zip
ZIP_DST=$input_dir/wav.zip
OUT_DIR=$input_dir/wav_full
wget $DATA_LINK -O $ZIP_DST
unzip $ZIP_DST -d $OUT_DIR
rm $ZIP_DST

DATA_LINK=https://atlas-group.cs.brown.edu/data/full/jpg.zip
ZIP_DST=$input_dir/jpg_full.zip
OUT_DIR=$input_dir/jpg_full
wget $DATA_LINK -O $ZIP_DST
unzip $ZIP_DST -d $OUT_DIR
rm $ZIP_DST

DATA_LINK=https://atlas-group.cs.brown.edu/data/small/jpg.zip
ZIP_DST=$input_dir/jpg_small.zip
OUT_DIR=$input_dir/jpg_small
wget $DATA_LINK -O $ZIP_DST
unzip $ZIP_DST -d $OUT_DIR
rm $ZIP_DST
