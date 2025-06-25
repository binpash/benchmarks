#!/bin/bash

cd "$(dirname "$0")" || exit 1
URL='https://atlas.cs.brown.edu/data'
in="./inputs"

mkdir -p "$in"

curl "$URL/web-index/stopwords.txt" > "$in/stopwords.txt" 

# TODO: Transform dataset into a format that can be used offline
