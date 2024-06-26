#!/bin/bash

#set -e

wiki_archive="https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z"
BENCH_TOP=${BENCH_TOP:-$(git rev-parse --show-toplevel)}

# . "$BENCH_TOP/scripts/utils.sh"
sudo apt-get install unzip

[ "$1" = "-c" ] && rm -rf en/ *.7z *.tar  500.txt 1000.txt full small

setup_dataset() {
  rm -rf ../1-grams.txt ../2-grams.txt 
  
  ## Downloading the dataset needs to happen for both small and large
  if [[ ! -d ./en ]]; then
    # wget $wiki_archive || eexit "cannot fetch wikipedia"
    # 7za x wikipedia-en-html.tar.7z
    tar -xvf wikipedia-en-html.tar
    wget http://ndr.md/data/wikipedia/index.txt # || eexit "cannot fetch wiki indices"
    # It is actually OK if we don't have this index since we download the 500/1000 below
  fi

  if [ "$1" = "--small" ]; then
    # 500 entries
    wget http://pac-n4.csail.mit.edu:81/pash_data/small/web-index.small.zip
    unzip web-index.small.zip
    mv small/500.txt .
    rm -rf small web-index.small.zip
  elif [ "$1" = "--full" ]; then
    the default full
    1000 entries
    wget http://pac-n4.csail.mit.edu:81/pash_data/full/web-index.full.zip
    unzip web-index.full.zip
    mv full/1000.txt .
    rm -rf full web-index.full.zip
  fi
}

setup_dataset $1