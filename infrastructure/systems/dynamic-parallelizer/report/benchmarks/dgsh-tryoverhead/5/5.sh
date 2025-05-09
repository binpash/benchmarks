#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file example/commit-stats.sh
#
# SYNOPSIS Highlight misspelled words
# DESCRIPTION
# Highlight the words that are misspelled in the command's standard
# input.
# Demonstrates stream buffering, the avoidance of pass-through
# constructs to avoid deadlocks, and the use of named streams.
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

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)

cat $INPUT_FILE >"$file1"

# Find errors

# Obtain list of words in text
cat "$file1" | 
tr '[:upper:]' '[:lower:]' | 
sed 's/[^a-z]/\n/g' | 
grep -v '^$' | 
sort | 
uniq | 
grep -v '^$' > "$file2"

# Ensure dictionary is compatibly sorted
sort /usr/share/dict/words > "$file3"

# List errors as a set difference
comm -23 "$file2" "$file3" 
