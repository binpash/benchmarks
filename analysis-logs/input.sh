#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
DIR=$REPO_TOP/analysis-logs/input
mkdir -p $DIR

# Set up Kaggle API
if [[ ! -d ~/.kaggle ]]; then
    mkdir ~/.kaggle
    echo "Place your kaggle.json in the ~/.kaggle directory."
fi
chmod 600 ~/.kaggle/kaggle.json

cd $DIR
kaggle datasets download -d eliasdabbas/web-server-access-logs
unzip web-server-access-logs
rm -f web-server-access-logs.zip client_hostname.csv
