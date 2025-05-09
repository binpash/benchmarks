#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file example/commit-stats.sh
#
# SYNOPSIS Git commit statistics
# DESCRIPTION
# Process the git history, and list the authors and days of the week
# ordered by the number of their commits.
# Demonstrates streams and piping through a function.
#
# Adjusted to not contain functions
#
# Recommended inputs:
# https://github.com/torvalds/linux
# https://github.com/kelsny/overcommitted
# https://github.com/996icu/996.ICU
#
#  Copyright 2012-2013 Diomidis Spinellis
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
REPO_DIR="$REPO_DIR"
OUTPUT_DIR="$OUTPUT_DIR"

file1="$OUTPUT_DIR/file1.txt"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i git -C $REPO_DIR log --format="%an:%ad" --date=default >"$file1"
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Authors ordered by number of commits"
# Order by frequency
logfile=$(generate_unique_file)
cmd='{print $1}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F: "$cmd" <"$file1" | sort | uniq | sort -rn

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Days ordered by number of commits"
# Order by frequency
logfile=$(generate_unique_file)
cmd='{print substr($2, 1, 3)}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F: "$cmd" <"$file1"  | sort | uniq | sort -rn
