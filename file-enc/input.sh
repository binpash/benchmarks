#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
input_dir="${eval_dir}/input"

DATA_LINK=https://atlas-group.cs.brown.edu/data/pcaps.zip
ZIP_DST=$input_dir/pcaps.zip
wget --no-check-certificate $DATA_LINK -O $ZIP_DST
unzip $ZIP_DST -d $input_dir
rm $ZIP_DST
