#!/bin/bash

OUTPUT_DIR="$OUTPUT_DIR"
extract_text="$SCRIPT_DIR/extract_text.sh"
bigrams_aux="$SCRIPT_DIR/bigrams_aux.sh"
trigrams_aux="$SCRIPT_DIR/trigrams_aux.sh"

cat $INPUT_FILE |
  sed "s#^#$WIKI/#" |
  $extract_text |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
  $WEB_INDEX_DIR/stem-words.js |
  tee "$OUTPUT_DIR/3grams" "$OUTPUT_DIR/2grams" "$OUTPUT_DIR/1grams"

cat "$OUTPUT_DIR/1grams" |
    sort |
    uniq -c |
    sort -rn > "$OUTPUT_DIR/1-grams.txt"

cat "$OUTPUT_DIR/2grams" |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    $bigrams_aux |
    sort |
    uniq -c |
    sort -rn > "$OUTPUT_DIR/2-grams.txt"

cat "$OUTPUT_DIR/3grams" |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    $trigrams_aux |
    sort |
    uniq -c |
    sort -rn > "$OUTPUT_DIR/3-grams.txt"
