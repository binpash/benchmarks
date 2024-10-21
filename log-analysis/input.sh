#!/bin/bash

set -e

# creates input/pcaps and input/nginx-logs

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/input"

DATA_LINK=https://atlas-group.cs.brown.edu/data/pcaps.zip
ZIP_DST=$input_dir/pcaps.zip
wget --no-check-certificate $DATA_LINK -O $ZIP_DST
unzip $ZIP_DST -d $input_dir
rm $ZIP_DST

DATA_LINK=https://atlas-group.cs.brown.edu/data/nginx.zip
ZIP_DST=$input_dir/nginx.zip
wget --no-check-certificate $DATA_LINK -O $ZIP_DST
unzip $ZIP_DST -d $input_dir
rm $ZIP_DST
