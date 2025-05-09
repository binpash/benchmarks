#!/bin/bash

# Adapted from: dspinellis/dgsh
# Source file example/web-log-report.sh
#
# SYNOPSIS Web log reporting
# DESCRIPTION
# Creates a report for a fixed-size web log file read from the standard input.
# Demonstrates the combined use of stores and named streams,
# the use of shell group commands and functions in the scatter block, and
# the use of cat(1) as a way to sequentially combine multiple streams.
# Used to measure throughput increase achieved through parallelism.
#
#
# Recommended inputs:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/3QBYB5
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

# Consitent sorting across machines

# Consistent sorting
export LC_ALL=C

export TZ=$(date +%Z)

# Print initial header only if DGSH_DRAW_EXIT is not set
if [ -z "${DGSH_DRAW_EXIT}" ]
then
    cat <<EOF
			WWW server statistics
			=====================

Summary
=======
EOF
fi

## Initialize temporary files
file_initial="$OUTPUT_DIR/file1.txt"
file_bytes="$OUTPUT_DIR/file2.txt"
file_hosts="$OUTPUT_DIR/file3.txt"
file_sorted_hosts="$OUTPUT_DIR/file4.txt"
file_unique_hosts="$OUTPUT_DIR/file5.txt"
file_domains="$OUTPUT_DIR/file6.txt"
file_requests="$OUTPUT_DIR/file7.txt"
file_times="$OUTPUT_DIR/file8.txt"
file_day_count="$OUTPUT_DIR/file9.txt"
file_dates="$OUTPUT_DIR/file10.txt"

# This file will capture a large portion of the processed data to be reused in subsequent parts
cat $INPUT_FILE > "$file_initial"

# Number of accesses
echo -n 'Number of accesses: '
wc -l < "$file_initial"

# Total transferred bytes
awk '{s += $NF} END {print s}' "$file_initial" > "$file_bytes"
echo -n 'Number of Gbytes transferred: '
awk '{print $1 / 1024 / 1024 / 1024}' "$file_bytes"

# Process Host names
awk '{print $1}' "$file_initial" > "$file_hosts"

# Number of accesses
echo -n 'Number of accesses: '
wc -l < "$file_hosts"

# Sorted hosts
sort "$file_hosts" > "$file_sorted_hosts"

# Unique hosts
uniq "$file_sorted_hosts" > "$file_unique_hosts"
echo -n 'Number of hosts: '
wc -l < "$file_unique_hosts"

# Number of TLDs
awk -F. '$NF !~ /[0-9]/ {print $NF}' "$file_unique_hosts" | sort -u | wc -l
echo -n 'Number of top level domains: '

# Top 10 hosts
echo
echo "Top 10 Hosts"
echo "Top 10 Hosts" | sed 's/./-/g'

uniq -c "$file_sorted_hosts" | sort -rn | head -10
echo

# Top 20 TLDs
echo
echo "Top 20 Level Domain Accesses"
echo "Top 20 Level Domain Accesses" | sed 's/./-/g'

awk -F. '$NF !~ /^[0-9]/ {print $NF}' "$file_sorted_hosts" | sort | uniq -c | sort -rn | head -20
echo

# Domains
awk -F. 'BEGIN {OFS = "."} $NF !~ /^[0-9]/ {$1 = ""; print}' "$file_sorted_hosts" | sort > "$file_domains"

# Number of domains
echo -n 'Number of domains: '
uniq "$file_domains" | wc -l

# Top 10 domains
echo
echo "Top 10 domains"
echo "Top 10 domains" | sed 's/./-/g'
uniq -c "$file_domains" | sort -rn | head -10 < "$file_domains"

# Hosts by volume
echo
echo "Top 10 Hosts by Transfer"
echo "Top 10 Hosts by Transfer" | sed 's/./-/g'
awk '    {bytes[$1] += $NF}
END {for (h in bytes) print bytes[h], h}' "$file_initial" | sort -rn | head -10

# Sorted page name requests
awk '{print $7}' "$file_initial" | sort > "$file_requests"

# Top 20 area requests (input is already sorted)
echo
echo "Top 20 area requests"
echo "Top 20 area requests" | sed 's/./-/g'
awk -F/ '{print $2}' "$file_requests" | uniq -c | sort -rn | head -20
# Number of different pages
echo -n 'Number of different pages: '
uniq "$file_requests" | wc -l

# Top 20 requests
echo
echo "Top 20 requests"
echo "Top 20 requests" | sed 's/./-/g'
uniq -c "$file_requests" | sort -rn | head -20

# Access time: dd/mmm/yyyy:hh:mm:ss
awk '{print substr($4, 2)}' "$file_initial" > "$file_times"

# Just dates
awk -F: '{print $1}' "$file_times" > "$file_dates"

# Number of days
echo -n 'Accesses per day: '
uniq "$file_dates" | wc -l > "$file_day_count"
awk '
BEGIN {
    getline NACCESS < "'"$file_initial"'"
}
{print NACCESS / $1}' "$file_day_count"

echo -n 'MBytes per day: '
awk '
BEGIN {
    getline NXBYTES < "'"$file_bytes"'"
}
{print NXBYTES / $1 / 1024 / 1024}' "$file_day_count"

echo
echo "Accesses by Date"
echo "Accesses by Date" | sed 's/./-/g'
uniq -c < "$file_dates"

# Accesses by day of week
echo
echo "Accesses by Day of Week"
echo "Accesses by Day of Week" | sed 's/./-/g'
sed 's|/|-|g' "$file_dates" | date -f - +%a 2>/dev/null | sort | uniq -c | sort -rn

# Accesses by Local Hour
echo
echo "Accesses by Local Hour"
echo "Accesses by Local Hour" | sed 's/./-/g'
awk -F: '{print $2}' "$file_times" | sort | uniq -c
