#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file example/reorder-columns.sh
# SYNOPSIS Reorder columns
# DESCRIPTION
# Reorder columns in a CSV document.
# Demonstrates the combined use of tee, cut, and paste.
#
#
# Recommended inputs:
# 5GB:   https://testfile.org/files-5GB
# 1GB:   https://bit.ly/1GB-testfile
# 700MB: http://aiweb.cs.washington.edu/research/projects/xmltk/xmldata/data/pir/psd7003.xml
# 120MB: http://aiweb.cs.washington.edu/research/projects/xmltk/xmldata/data/dblp/dblp.xml
# 23MB:  http://aiweb.cs.washington.edu/research/projects/xmltk/xmldata/data/nasa/nasa.xml
# 2MB:   http://aiweb.cs.washington.edu/research/projects/xmltk/xmldata/data/courses/uwm.xml
#
#
#  Copyright 2016 Marios Fragkoulis
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

# Extract columns 5 and 6, save to temp1
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cut -d ',' -f 5-6 "$INPUT_FILE" > "$file1"

# Extract columns 2, 3, and 4, save to temp2
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cut -d ',' -f 2-4 "$INPUT_FILE" > "$file2"

# Combine the columns
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i paste -d ',' "$file1" "$file2"
