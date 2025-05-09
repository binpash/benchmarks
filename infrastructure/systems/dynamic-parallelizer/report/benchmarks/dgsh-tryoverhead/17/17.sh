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

# Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)
file4=$(mktemp)

cat $INPUT_FILE > $file1

# Extract columns 5 and 6, save to temp1
cut -d ',' -f 5-6 "$file1" > "$file2"

# Extract columns 2, 3, and 4, save to temp2
cut -d ',' -f 2-4 "$file1" > "$file3"

# Combine the columns
paste -d ',' "$file2" "$file3"
