#!/bin/bash

cat $INPUT_FILE |
sed "s#^#$WEB_INDEX_DIR#" |
iconv -c -t ascii//TRANSLIT |
pandoc +RTS -K64m -RTS --from html --to plain --quiet |
tr -cs A-Za-z '\n' |
tr A-Z a-z |
grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
$WEB_INDEX_DIR/stem-words.js
