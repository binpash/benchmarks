REPO_TOP=$(git rev-parse --show-toplevel)
export TEST_BASE=$REPO_TOP/web-index
export SCRIPT_DIR="$TEST_BASE"/scripts
export WEB_INDEX_DIR="$TEST_BASE"/inputs
export WIKI="$TEST_BASE"/inputs/articles

cd $(dirname "$0") || exit 1

output_base="$1"

rm -f {1,2,3}grams
mkfifo {1,2,3}grams

extract_text="$SCRIPT_DIR/extract_text.sh"
bigrams_aux="$SCRIPT_DIR/bigrams_aux.sh"
trigrams_aux="$SCRIPT_DIR/trigrams_aux.sh"

cat "$INPUT_FILE" |
  sed "s#^#$WIKI/#" |
  $extract_text |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf "$WEB_INDEX_DIR/stopwords.txt" |
  "$SCRIPT_DIR/stem-words.js" |
  tee 3grams 2grams 1grams > /dev/null &

cat 1grams |
    sort |
    uniq -c |
    sort -rn > "$output_base/1-grams.txt" &

cat 2grams |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    $bigrams_aux |
    sort |
    uniq -c |
    sort -rn > "$output_base/2-grams.txt" &

cat 3grams |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    $trigrams_aux |
    sort |
    uniq -c |
    sort -rn > "$output_base/3-grams.txt"

rm -f {1,2,3}grams
