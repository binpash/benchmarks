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

## Initialize the necessary temporary files
file1=$(mktemp)

git -C $REPO_DIR log --format="%an:%ad" --date=default >"$file1"
echo "Authors ordered by number of commits"
# Order by frequency
awk -F: '{print $1}' <"$file1" | sort | uniq | sort -rn

echo "Days ordered by number of commits"
# Order by frequency
awk -F: '{print substr($2, 1, 3)}' <"$file1"  | sort | uniq | sort -rn

