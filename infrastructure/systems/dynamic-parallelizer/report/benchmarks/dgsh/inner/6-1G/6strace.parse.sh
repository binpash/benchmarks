#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file example/word-properties.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Word properties
# DESCRIPTION
# Read text from the standard input and list words
# containing a two-letter palindrome, words containing
# four consonants, and words longer than 12 characters.
#
# Demonstrates the use of paste as a gather function
#
# Recommended inputs:
# ftp://sunsite.informatik.rwth-aachen.de/pub/mirror/ibiblio/gutenberg/1/3/139/139.txt
# https://ocw.mit.edu/ans7870/6/6.006/s08/lecturenotes/files/t8.shakespeare.txt
#
#  Copyright 2013 Diomidis Spinellis
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
generate_unique_file() {
    local dir="$OUTPUT_DIR"
    local prefix="strace_log"
    local counter_file="$dir/${prefix}"
    if [ ! -f "$counter_file" ]; then
        echo 0 > "$counter_file"
    fi
    local counter
    counter=$(<"$counter_file")
    counter=$((counter + 1))
    echo "$counter" > "$counter_file"
    local filename="$dir/${prefix}_${counter}"
    echo "$filename"
}

export STRACE="strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i"

INPUT_FILE="$INPUT_FILE"
OUTPUT_DIR="$OUTPUT_DIR"

# Output files
file1="$OUTPUT_DIR/file1.txt"
file2="$OUTPUT_DIR/file2.txt"
file3="$OUTPUT_DIR/file3.txt"
file4="$OUTPUT_DIR/file4.txt"
file5="$OUTPUT_DIR/file5.txt"

# Stream input from file and split input one word per line
# Create list of unique words
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i tr -cs a-zA-Z '\n' < "$INPUT_FILE" | sort -u > "$file1"
$PARSE $logfile > $(generate_unique_file)

# List two-letter palindromes
logfile=$(generate_unique_file)
expr='s/.*\(.\)\(.\)\2\1.*/p: \1\2-\2\1/;t;g'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i sed "$expr" "$file1" > "$file2"
$PARSE $logfile > $(generate_unique_file)

# List four consecutive consonants
logfile=$(generate_unique_file)
expr='s/.*([^aeiouyAEIOUY]{4}).*/c: \1/;t;g'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i sed -E "$expr" "$file1" > "$file3"
$PARSE $logfile > $(generate_unique_file)

# List length of words longer than 12 characters
logfile=$(generate_unique_file)
cmd='{if (length($1) > 12) print "l:", length($1);
	else print ""}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk $cmd "$file1" > "$file4"
$PARSE $logfile > $(generate_unique_file)

# Paste the four streams side-by-side
# List only words satisfying one or more properties
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i paste "$file1" "$file2" "$file3" "$file4" | fgrep :
$PARSE $logfile > $(generate_unique_file)
