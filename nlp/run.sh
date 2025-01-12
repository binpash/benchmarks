#!/bin/bash --posix

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd "$SUITE_DIR"

if [[ "$1" == "--small" ]]; then
    export ENTRIES=10
    export IN="$SUITE_DIR/inputs/pg-small"
else
    export ENTRIES=120
    export IN="$SUITE_DIR/inputs/pg"
fi

export BENCHMARK_CATEGORY="nlp"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

mkdir -p "outputs"

# Define the script names in a single variable
script_names="syllable_words_1
syllable_words_2
letter_words
bigrams_appear_twice
bigrams
compare_exodus_genesis
count_consonant_seq
count_morphs
count_trigrams
count_vowel_seq
count_words
find_anagrams
merge_upper
sort
sort_words_by_folding
sort_words_by_num_of_syllables
sort_words_by_rhyming
trigram_rec
uppercase_by_token
uppercase_by_type
verses_2om_3om_2instances
vowel_sequencies_gr_1K
words_no_vowels"

mkdir -p "outputs"

echo "Executing nlp $(date)"

# Loop through each script name from the variable
while IFS= read -r script; do
    script_file="./scripts/$script.sh"
    output_dir="./outputs/$script/"

    mkdir -p "$output_dir"

    export BENCHMARK_SCRIPT="$(realpath "$script_file")"
    echo "$script"
    $BENCHMARK_SHELL "$script_file" "$output_dir"
    echo "$?"
done <<< "$script_names"
