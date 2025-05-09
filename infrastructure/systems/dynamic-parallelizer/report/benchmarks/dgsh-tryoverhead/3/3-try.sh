#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file: example/code-metrics.sh
#
# SYNOPSIS C code metrics
# DESCRIPTION
# Process a directory containing C source code, and produce a summary
# of various metrics.
# Demonstrates nesting, commands without input.
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

try=/srv/hs/deps/try/try

## Initialize the necessary temporary files
$try -y file1=$(mktemp)
$try -y file2=$(mktemp)
$try -y file3=$(mktemp)
$try -y file4=$(mktemp)

$try -y find "$REPO_DIR" \( -name \*.c -or -name \*.h \) -type f -print0 >"$file1"

$try -y echo -n 'FNAMELEN: '

$try -y "tr \\0 \\n <"$file1" |
sed 's|^.*/||' |
awk '{s += length($1); n++} END {
    if (n>0)
        print s / n;
    else
        print 0; }'"

$try -y xargs -0 /bin/cat <"$file1" >"$file2"

sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" <"$file2" |
    cpp -P >"$file3"

# Structure definitions
echo -n 'NSTRUCT: '

egrep -c 'struct[   ]*{|struct[   ]*[a-zA-Z_][a-zA-Z0-9_]*[       ]*{' <"$file3"
#}} (match preceding openings)

# Type definitions
echo -n 'NTYPEDEF: '
grep -cw typedef <"$file3"

# Use of void
echo -n 'NVOID: '
grep -cw void <"$file3"

# Use of gets
echo -n 'NGETS: '
grep -cw gets <"$file3"

# Average identifier length
echo -n 'IDLEN: '

tr -cs 'A-Za-z0-9_' '\n' <"$file3" |
sort -u |
awk '/^[A-Za-z]/ { len += length($1); n++ } END {
    if (n>0)
        print len / n;
    else
        print 0; }'

echo -n 'CHLINESCHAR: '
wc -lc  <"$file2" |
    awk '{OFS=":"; print $1, $2}'

echo -n 'NCCHAR: '
sed 's/#/@/g' <"$file2" |
cpp -traditional -P |
wc -c |
awk '{OFMT = "%.0f"; print $1/1000}'

# Number of comments
echo -n 'NCOMMENT: '
egrep -c '/\*|//' <"$file2"

# Occurences of the word Copyright
echo -n 'NCOPYRIGHT: '
grep -ci copyright <"$file2"

# C files
find "$@" -name \*.c -type f -print0 >"$file2"

# Convert to newline separation for counting
tr \\0 \\n <"$file2" >"$file3"

# Number of C files
echo -n 'NCFILE: '
wc -l <"$file3"

# Number of directories containing C files
echo -n 'NCDIR: '
sed 's,/[^/]*$,,;s,^.*/,,' <"$file3" |
sort -u |
wc -l

# C code
xargs -0 /bin/cat <"$file2" >"$file3"

# Lines and characters
echo -n 'CLINESCHAR: '
wc -lc <"$file3" |
awk '{OFS=":"; print $1, $2}'

# C code without comments and strings
sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" <"$file3" |
cpp -P >"$file4"

# Number of functions
echo -n 'NFUNCTION: '
grep -c '^{' <"$file4"

# Number of gotos
echo -n 'NGOTO: '
grep -cw goto <"$file4"

# Occurrences of the register keyword
echo -n 'NREGISTER: '
grep -cw register <"$file4"

# Number of macro definitions
echo -n 'NMACRO: '
grep -c '@[   ]*define[   ][   ]*[a-zA-Z_][a-zA-Z0-9_]*(' <"$file4"
# Number of include directives
echo -n 'NINCLUDE: '
grep -c '@[   ]*include' <"$file4"

# Number of constants
echo -n 'NCONST: '
grep -ohw '[0-9][x0-9][0-9a-f]*' <"$file4" | wc -l 


# Header files
echo -n 'NHFILE: '
find "$@" -name \*.h -type f |
wc -l
