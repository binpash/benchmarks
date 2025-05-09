#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file: example/compress-compare.sh
#
# SYNOPSIS Compression benchmark
# DESCRIPTION
# Report file type, length, and compression performance for
# data received from the standard input.  The data never touches the
# disk.
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
echo $INPUT_FILE
file1=$(mktemp)
cat $INPUT_FILE >"$file1"
printf 'File type:\t'
file - <"$file1"

printf 'Original size:\t'
wc -c <"$file1"

printf 'xz:\t\t'
xz -c <"$file1" | wc -c

printf 'bzip2:\t\t'
bzip2 -c <"$file1" | wc -c

printf 'gzip:\t\t'
gzip -c <"$file1" | wc -c
