#!/bin/sh

cd "$(dirname "$0")"

directory_path="articles"

if [ ! -d "$directory_path" ]; then
	    echo "Error: Directory does not exist."
	        exit 1
fi

# ensure a local ./tmp directory exists for sorting
mkdir -p ./tmp
export TMPDIR=./tmp

# find all files, remove prefix, sort them, and write to a text file
find "$directory_path" -type f | sed 's|./wikipedia/en/articles/||' | sort > index.txt

echo "File paths have been saved to all_files_paths.txt"


# default values
WINDOW=5
TARGET="sh-only"
LOG="enable"
INPUT="small"

# parse arguments
while [ "$#" -gt 0 ]; do
    case $1 in
        --window) WINDOW="$2"; shift ;;
        --target) TARGET="$2"; shift ;;
        --log) LOG="$2"; shift ;;
        --input) INPUT="$2"; shift ;;
    esac
    shift
done

# setz environment variables
TEST_BASE=$(dirname "$(realpath "$0")")
echo "TEST_BASE: $TEST_BASE"
LOCAL_NAME=$(basename "$TEST_BASE")
echo "LOCAL_NAME: $LOCAL_NAME"
OUTPUT_BASE="$TEST_BASE/output/$LOCAL_NAME"

if [ "$INPUT" = "small" ]; then
    export INPUT_FILE="$TEST_BASE/index_small.txt"
elif [ "$INPUT" = "full" ]; then
    export INPUT_FILE="$TEST_BASE/index.txt"
fi
export WEB_INDEX_DIR="$TEST_BASE"
export SCRIPT_DIR="$TEST_BASE"
export WIKI="$TEST_BASE/articles"

do_cleanup() {
    echo "Cleaning up..."
    rm -f *grams *-grams.txt
}

do_run() {
    output_base=$1
    start_time=$(date +%s)
    
    echo "Running integrated script"
    
    mkfifo {1,2,3}grams
    
    extract_text="$SCRIPT_DIR/extract_text.sh"
    bigrams_aux="$SCRIPT_DIR/bigrams_aux.sh"
    trigrams_aux="$SCRIPT_DIR/trigrams_aux.sh"
    
    cat "$INPUT_FILE" |
      sed "s#^#$WIKI/#" | head |
      $extract_text |
      tr -cs A-Za-z '\n' |
      tr A-Z a-z |
      grep -vwFf "$WEB_INDEX_DIR/stopwords.txt" |
      "$WEB_INDEX_DIR/stem-words.js" |
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
    
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "$duration" > "$output_base/sh_time"
}

# create output directory
mkdir -p "$OUTPUT_BASE"

# run the integrated process
do_cleanup
do_run "$OUTPUT_BASE"
