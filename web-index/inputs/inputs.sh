#!/bin/bash

#set -e

wiki_archive="https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z"


rm-files(){
  echo "${@}"
  rm -r "${@}"
  exit 0
}

if [[ "$1" = "-c" ]]; then
    rm-files ../data/ *.7z *.tar  500.txt 1000.txt full small
fi

setup_dataset_web_index() {
#  rm -rf ../1-grams.txt ../2-grams.txt 
  
  ## Downloading the dataset needs to happen for both small and large
    if [[ ! -d ../data ]]; then
        mkdir ../data
        wget -P ../data $wiki_archive || eexit "cannot fetch wikipedia"
        7za x ../data/wikipedia-en-html.tar.7z
        tar -xvf ../data/wikipedia-en-html.tar
        wget -P ../data atlas-group.cs.brown.edu/data/wikipedia/index.txt 
    
    fi
} 

setup_dataset_web_index