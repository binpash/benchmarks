#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
export TEST_BASE=$REPO_TOP/web-index
export SCRIPT_DIR="$TEST_BASE/scripts"
export WEB_INDEX_DIR="$TEST_BASE/inputs"
export WIKI="$TEST_BASE/inputs/articles"

cd "$(dirname "$0")" || exit 1

output_base="$1"

extract_text="$SCRIPT_DIR/extract_text.sh"
bigrams_aux="$SCRIPT_DIR/bigrams_aux.sh"
trigrams_aux="$SCRIPT_DIR/trigrams_aux.sh"

# Process the input file and generate n-grams
sed "s#^#$WIKI/#" "$INPUT_FILE" |
  "$extract_text" |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf "$WEB_INDEX_DIR/stopwords.txt" |
  "$SCRIPT_DIR/stem-words.js" |
  tee >(sort | uniq -c | sort -rn > "$output_base/1-grams.txt") \
      >(tr -cs A-Za-z '\n' | tr A-Z a-z | "$bigrams_aux" | sort | uniq -c | sort -rn > "$output_base/2-grams.txt") \
      >(tr -cs A-Za-z '\n' | tr A-Z a-z | "$trigrams_aux" | sort | uniq -c | sort -rn > "$output_base/3-grams.txt") > /dev/null