#!/bin/bash

OUTPUT_DIR="$OUTPUT_DIR"
extract_text="$SCRIPT_DIR/extract_text.sh"
bigrams_aux="$SCRIPT_DIR/bigrams_aux.sh"
trigrams_aux="$SCRIPT_DIR/trigrams_aux.sh"

bigrams_aux()
{
    ( mkfifo s2 > /dev/null ) ;
    ( mkfifo s3 > /dev/null ) ;

    sed '$d' s2 > s3 &
    tee s2 |
        tail +2 |
        paste s3 -
    rm s2
    rm s3
}

trigrams_aux()
{
    s2=$(mktemp -u)
    s3=$(mktemp -u)

    mkfifo $s2 $s3

    tee $s2 |
        tail +2 |
        paste $s2 - |
        tee $s3 |
        cut -f 1 |
        tail +3 |
        paste $s3 - |
        sed "\$d" |
        sed "\$d"

    rm $s2 $s3
}


extract_text()
{
    while read -r line
    do
        cat $line |
            iconv -c -t ascii//TRANSLIT |
            pandoc +RTS -K64m -RTS --from html --to plain --quiet
    done
}

cat $INPUT_FILE | 
  sed "s#^#$WIKI/#" |
  extract_text |
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
    bigrams_aux |
    sort |
    uniq -c |
    sort -rn > "$OUTPUT_DIR/2-grams.txt"

cat "$OUTPUT_DIR/3grams" |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    trigrams_aux |
    sort |
    uniq -c |
    sort -rn > "$OUTPUT_DIR/3-grams.txt"
