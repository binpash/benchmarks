#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file example/text-properties.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Text properties
# DESCRIPTION
# Read text from the standard input and create files
# containing word, character, digram, and trigram frequencies.
#
# Demonstrates the use of scatter blocks without output and the use
# of stores within the scatter block.
#
# Functions have been removed.
#
# Example:
# curl ftp://sunsite.informatik.rwth-aachen.de/pub/mirror/ibiblio/gutenberg/1/3/139/139.txt | text-properties
#
#  Copyright 2013 Diomidis Spinellis
#
# Recommended inputs:
# ftp://sunsite.informatik.rwth-aachen.de/pub/mirror/ibiblio/gutenberg/1/3/139/139.txt
# https://ocw.mit.edu/ans7870/6/6.006/s08/lecturenotes/files/t8.shakespeare.txt
#
#
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

# Input and Output directories from environment variables
INPUT_FILE="$INPUT_FILE"
OUTPUT_DIR="$OUTPUT_DIR"

# Output files
file1="$OUTPUT_DIR/file1.txt"
file2="$OUTPUT_DIR/file2.txt"

# Store number of characters to use in awk below
nchars=$(wc -c < "$INPUT_FILE")


# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Split input one word per line
tr -cs a-zA-Z '\n' < "$INPUT_FILE" > "$file1"

# Digram frequency
echo "Digram frequency"
perl -ne 'for ($i = 0; $i < length($_) - 2; $i++) {
	print substr($_, $i, 2), "\n";
}' < "$file1" |
awk '{count[$1]++} END {for (i in count) print count[i], i}' |
sort -rn

# Trigram frequency
echo "Trigram frequency"
perl -ne 'for ($i = 0; $i < length($_) - 3; $i++) {
	print substr($_, $i, 3), "\n";
}' < "$file1" |
awk '{count[$1]++} END {for (i in count) print count[i], i}' |
sort -rn

# Word frequency
echo "Word frequency"
awk '{count[$1]++} END {for (i in count) print count[i], i}' < "$file1" |
sort -rn

# Character frequency
# Print absolute
echo "Character frequency"
sed 's/./&\n/g' < "$INPUT_FILE" |
awk '{count[$1]++} END {for (i in count) print count[i], i}' |
sort -rn | tee "$file2"

# Print relative
echo "Relative character frequency"
awk -v NCHARS="$nchars" 'BEGIN {
		OFMT = "%.2g%%"}
		{print $1, $2, $1 / NCHARS * 100}' "$file2"
