#!/bin/bash
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
WIKI=${WIKI:-$PASH_TOP/web-index}

cat $WIKI/input/index.txt |
sed "s#^#$WIKI#" |
iconv -c -t ascii//TRANSLIT |
pandoc +RTS -K64m -RTS --from html --to plain --quiet |
tr -cs A-Za-z '\n' |
tr A-Z a-z |
grep -vwFf $WIKI/stopwords.txt |
$WIKI/stem-words.js

