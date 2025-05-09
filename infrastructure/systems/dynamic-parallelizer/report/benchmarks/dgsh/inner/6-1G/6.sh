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
tr -cs a-zA-Z '\n' < "$INPUT_FILE" |
sort -u > "$file1"

# List two-letter palindromes
sed 's/.*\(.\)\(.\)\2\1.*/p: \1\2-\2\1/;t;g' "$file1" > "$file2"

# List four consecutive consonants
sed -E 's/.*([^aeiouyAEIOUY]{4}).*/c: \1/;t;g' "$file1" > "$file3"

# List length of words longer than 12 characters
awk '{if (length($1) > 12) print "l:", length($1);
	else print ""}' "$file1" > "$file4"

# Paste the four streams side-by-side
# List only words satisfying one or more properties
paste "$file1" "$file2" "$file3" "$file4" | 
fgrep :
