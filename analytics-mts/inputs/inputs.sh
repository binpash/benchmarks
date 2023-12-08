#!/bin/bash

rm-files(){
  echo "${@}"
  rm -r "${@}"
  exit 0
}

if [[ "$1" = "-c" ]]; then
    rm-files ../data/ *.bz2 'in.csv' 'in_small.csv'
    exit
fi


setup_dataset_mts() {
    if [[ ! -d ../data ]]; then
        mkdir -p ../data
    fi
    if [ ! -f ../data/in.csv ] && [ "$1" = "--full" ]; then
        # yesterday=$(date --date='1 days ago' +'%y-%m-%d')
        # curl httpsci://www.balab.aueb.gr/~dds/oasa-$yesterday.bz2 |
        curl 'https://www.balab.aueb.gr/~dds/oasa-2021-01-08.bz2' | bzip2 -d > ../data/in.csv
    if [ $? -ne 0 ]; then
        echo "oasa-2021-01-08.bz2 / bzip2 not available, contact the pash authors"
        #exit 1
    fi
    elif [ ! -f ../data/in_small.csv ] && [ "$1" = "--small" ]; then
        if [ ! -f ../data/in_small.csv ]; then                                                       
        echo "Generating small-size inputs"                                                  
        # FIXME PR: Do we need all of them?                                                  
        curl -sf 'http://pac-n4.csail.mit.edu:81/pash_data/small/in_small.csv' > ../data/in_small.csv
    fi                                                                                     
  fi
}

setup_dataset_mts $1