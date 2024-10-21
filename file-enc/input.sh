#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
input_dir="${eval_dir}/input"

mkdir -p $input_dir

DATA_LINK=https://atlas-group.cs.brown.edu/data/pcaps.zip
ZIP_DST=$input_dir/pcaps.zip

if [ ! -d $input_dir/pcaps ]; then
    wget --no-check-certificate $DATA_LINK -O $ZIP_DST
    unzip $ZIP_DST -d $input_dir
    rm $ZIP_DST
fi
