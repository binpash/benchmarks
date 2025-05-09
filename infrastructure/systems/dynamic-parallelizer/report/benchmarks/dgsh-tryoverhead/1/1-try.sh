#!/bin/bash
try=try

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
$try -y echo $INPUT_FILE
$try -y file1=$(mktemp)
$try -y cat $INPUT_FILE >"$file1"
$try -y printf 'File type:\t'
$try -y file - <"$file1"

$try -y printf 'Original size:\t'
$try -y wc -c <"$file1"

$try -y printf 'xz:\t\t'
$try -y xz -c <"$file1" | $try -y wc -c
$try -y "xz -c <\"$file1\" | wc -c"

$try -y printf 'bzip2:\t\t'
$try -y bzip2 -c <"$file1" | $try -y wc -c

$try -y printf 'gzip:\t\t'
$try -y gzip -c <"$file1" | $try -y wc -c
