#!/bin/bash
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
IN=${IN:-$PASH_TOP/web-index/input/index.txt}
WEB_INDEX_DIR=${WEB_INDEX_DIR:-$PASH_TOP/web-index/input}
WIKI=${WIKI:-$PASH_TOP/web-index/}

bigrams_aux() {
    paste <(sed '$d') <(tail +2)
}

bigram_aux_map() {
    IN=$1
    OUT=$2
    AUX_HEAD=$3
    AUX_TAIL=$4

    head -n 1 "$IN" > "$AUX_HEAD" &
    tail -n 1 "$IN" > "$AUX_TAIL" &
    paste <(sed '$d' "$IN") <(tail +2 "$IN") > "$OUT" &
    wait
}

bigram_aux_reduce() {
    IN1=$1
    AUX_HEAD1=$2
    AUX_TAIL1=$3
    IN2=$4
    AUX_HEAD2=$5
    AUX_TAIL2=$6
    OUT=$7
    AUX_HEAD_OUT=$8
    AUX_TAIL_OUT=$9

    cp "$AUX_HEAD1" "$AUX_HEAD_OUT" &
    cp "$AUX_TAIL2" "$AUX_TAIL_OUT" &
    {
        cat "$IN1"
        paste "$AUX_TAIL1" "$AUX_HEAD2"
        cat "$IN2"
    } > "$OUT" &
    wait
}

trigrams_aux() {
    paste <(paste <(sed '$d') <(tail +2)) <(tail +3) | sed '$d' | sed '$d'
}

extract_text() {
    while read -r line; do
        iconv -c -t ascii//TRANSLIT < "$line" |
            pandoc +RTS -K64m -RTS --from html --to plain --quiet
    done
}

export -f extract_text

# Process input and generate n-grams
head "$IN" |
    sed "s#^#$WIKI#" |
    extract_text |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    grep -vwFf "$WIKI/stopwords.txt" |
    "$WIKI/stem-words.js" |
    tee >(tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn > 1-grams.txt) \
        >(tr -cs A-Za-z '\n' | tr A-Z a-z | bigrams_aux | sort | uniq -c | sort -rn > 2-grams.txt) \
        >(tr -cs A-Za-z '\n' | tr A-Z a-z | trigrams_aux | sort | uniq -c | sort -rn > 3-grams.txt) > /dev/null
