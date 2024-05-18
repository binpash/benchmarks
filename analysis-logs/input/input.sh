#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
DIR=$REPO_TOP/analysis-logs/input
mkdir -p $DIR
cd $DIR

if [[ $1 == "--kaggle" ]]; then
    # Set up Kaggle API
    if [[ ! -d ~/.kaggle ]]; then
        mkdir ~/.kaggle
        echo "Place your kaggle.json in the ~/.kaggle directory."
    fi
    chmod 600 ~/.kaggle/kaggle.json

    if [[ ! -f nginx.zip ]]; then
        kaggle datasets download -d eliasdabbas/web-server-access-logs
        unzip web-server-access-logs
        rm -f web-server-access-logs.zip client_hostname.csv
    fi
else
    if [[ ! -f nginx.zip ]]; then
        # TODO: replace with omega URL
        # wget -O nginx.zip "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/3QBYB5/NXKB6J"
        # unzip web-server-access-logs
        # rm -f web-server-access-logs.zip
        echo "Not implemented yet."
        exit 1
    fi
fi
