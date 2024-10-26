#!/bin/bash
mkfifo {1,2,3}grams

extract_text="$SCRIPT_DIR/extract_text.sh"
bigrams_aux="$SCRIPT_DIR/bigrams_aux.sh"
trigrams_aux="$SCRIPT_DIR/trigrams_aux.sh"

cat $INPUT_FILE |
  sed "s#^#$WIKI/#" | head |
  $extract_text |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
  $WEB_INDEX_DIR/stem-words.js |
  tee 3grams 2grams 1grams

cat 1grams |
    sort |
    uniq -c |
    sort -rn > 1-grams.txt

cat 2grams |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    $bigrams_aux |
    sort |
    uniq -c |
    sort -rn > 2-grams.txt

cat 3grams |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    $trigrams_aux |
    sort |
    uniq -c |
    sort -rn > 3-grams.txt

rm -f {1,2,3}grams
