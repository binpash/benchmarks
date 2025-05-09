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

# Input and Output directories from environment variables
INPUT_FILE="$INPUT_FILE"
OUTPUT_DIR="$OUTPUT_DIR"

# Output files
file1="$OUTPUT_DIR/file1.txt"
file2="$OUTPUT_DIR/file2.txt"


logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i tr -cs a-zA-Z '\n' < "$INPUT_FILE" > "$file1"
$PARSE $logfile > $(generate_unique_file)

# Digram frequency
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Digram frequency"
$PARSE $logfile > $(generate_unique_file)

logfile=$(generate_unique_file)
perl_cmd='for ($i = 0; $i < length($_) - 2; $i++) { print substr($_, $i, 2), "\n"; }'
cmd='{count[$1]++} END {for (i in count) print count[i], i}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i perl -ne "$cmd" < "$file1" | awk "$cmd" | sort -rn
$PARSE $logfile > $(generate_unique_file)

# Trigram frequency
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Trigram frequency"
$PARSE $logfile > $(generate_unique_file)

logfile=$(generate_unique_file)
perl_cmd='for ($i = 0; $i < length($_) - 3; $i++) { print substr($_, $i, 3), "\n"; }'
cmd='{count[$1]++} END {for (i in count) print count[i], i}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i perl -ne "$perl_cmd" < "$file1" | awk "$cmd" | sort -rn
$PARSE $logfile > $(generate_unique_file)

# Word frequency
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Word frequency"
$PARSE $logfile > $(generate_unique_file)

logfile=$(generate_unique_file)
cmd='{count[$1]++} END {for (i in count) print count[i], i}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" < "$file1" | sort -rn
$PARSE $logfile > $(generate_unique_file)

# Store number of characters to use in awk below
nchars=$(wc -c < "$INPUT_FILE")

# Character frequency
# Print absolute
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Character frequency"
$PARSE $logfile > $(generate_unique_file)

logfile=$(generate_unique_file)
cmd='{count[$1]++} END {for (i in count) print count[i], i}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i sed 's/./&\n/g' < "$INPUT_FILE" | awk "$cmd" | sort -rn | tee "$file2"
$PARSE $logfile > $(generate_unique_file)

# Print relative

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Relative character frequency"
$PARSE $logfile > $(generate_unique_file)

logfile=$(generate_unique_file)
cmd='-v NCHARS="$nchars" BEGIN { OFMT = "%.2g%%"} {print $1, $2, $1 / NCHARS * 100}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -v "$cmd" "$file2"
$PARSE $logfile > $(generate_unique_file)
