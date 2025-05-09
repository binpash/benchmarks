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
export TRY_CMD="$TRY -y"
# Consistent sorting
export LC_ALL=C

# Print initial header only if DGSH_DRAW_EXIT is not set
if [ -z "${DGSH_DRAW_EXIT}" ]
then
$TRY_CMD cat <<EOF
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
$TRY_CMD cat $INPUT_FILE > "$file_initial"

# Number of accesses
$TRY_CMD echo -n 'Number of accesses: '
$TRY_CMD wc -l < "$file_initial"

# Total transferred bytes
$TRY_CMD awk '{s += $NF} END {print s}' "$file_initial" > "$file_bytes"
$TRY_CMD echo -n 'Number of Gbytes transferred: '
$TRY_CMD awk '{print $1 / 1024 / 1024 / 1024}' "$file_bytes"

# Process Host names
$TRY_CMD awk '{print $1}' "$file_initial" > "$file_hosts"

# Number of accesses
$TRY_CMD echo -n 'Number of accesses: '
$TRY_CMD wc -l < "$file_hosts"

# Sorted hosts
$TRY_CMD sort "$file_hosts" > "$file_sorted_hosts"

# Unique hosts
$TRY_CMD uniq "$file_sorted_hosts" > "$file_unique_hosts"
$TRY_CMD echo -n 'Number of hosts: '
$TRY_CMD wc -l < "$file_unique_hosts"

# Number of TLDs
$TRY_CMD awk -F. '$NF !~ /[0-9]/ {print $NF}' "$file_unique_hosts" | sort -u | wc -l
$TRY_CMD echo -n 'Number of top level domains: '

# Top 10 hosts
$TRY_CMD echo
$TRY_CMD echo "Top 10 Hosts"
$TRY_CMD echo "Top 10 Hosts" | sed 's/./-/g'

$TRY_CMD uniq -c "$file_sorted_hosts" | sort -rn | head -10
$TRY_CMD echo

# Top 20 TLDs
$TRY_CMD echo
$TRY_CMD echo "Top 20 Level Domain Accesses"
$TRY_CMD echo "Top 20 Level Domain Accesses" | sed 's/./-/g'

$TRY_CMD awk -F. '$NF !~ /^[0-9]/ {print $NF}' "$file_sorted_hosts" | sort | uniq -c | sort -rn | head -20
$TRY_CMD echo

# Domains
$TRY_CMD awk -F. 'BEGIN {OFS = "."} $NF !~ /^[0-9]/ {$1 = ""; print}' "$file_sorted_hosts" | sort > "$file_domains"

# Number of domains
$TRY_CMD echo -n 'Number of domains: '
$TRY_CMD uniq "$file_domains" | wc -l

# Top 10 domains
$TRY_CMD echo
$TRY_CMD echo "Top 10 domains"
$TRY_CMD echo "Top 10 domains" | sed 's/./-/g'
$TRY_CMD uniq -c "$file_domains" | sort -rn | head -10 < "$file_domains"

# Hosts by volume
$TRY_CMD echo
$TRY_CMD echo "Top 10 Hosts by Transfer"
$TRY_CMD echo "Top 10 Hosts by Transfer" | sed 's/./-/g'
$TRY_CMD awk '    {bytes[$1] += $NF}
END {for (h in bytes) print bytes[h], h}' "$file_initial" | sort -rn | head -10

# Sorted page name requests
$TRY_CMD awk '{print $7}' "$file_initial" | sort > "$file_requests"

# Top 20 area requests (input is already sorted)
$TRY_CMD echo
$TRY_CMD echo "Top 20 area requests"
$TRY_CMD echo "Top 20 area requests" | sed 's/./-/g'
$TRY_CMD awk -F/ '{print $2}' "$file_requests" | uniq -c | sort -rn | head -20
# Number of different pages
$TRY_CMD echo -n 'Number of different pages: '
$TRY_CMD uniq "$file_requests" | wc -l

# Top 20 requests
$TRY_CMD echo
$TRY_CMD echo "Top 20 requests"
$TRY_CMD echo "Top 20 requests" | sed 's/./-/g'
$TRY_CMD uniq -c "$file_requests" | sort -rn | head -20

# Access time: dd/mmm/yyyy:hh:mm:ss
$TRY_CMD awk '{print substr($4, 2)}' "$file_initial" > "$file_times"

# Just dates
$TRY_CMD awk -F: '{print $1}' "$file_times" > "$file_dates"

# Number of days
$TRY_CMD echo -n 'Accesses per day: '
$TRY_CMD uniq "$file_dates" | wc -l > "$file_day_count"
$TRY_CMD awk 'BEGIN {getline NACCESS < "'"$file_initial"'"}{print NACCESS / $1}' "$file_day_count"
$TRY_CMD echo -n 'MBytes per day: '
$TRY_CMD awk 'BEGIN {getline NXBYTES < "'"$file_bytes"'"}{print NXBYTES / $1 / 1024 / 1024}' "$file_day_count"

$TRY_CMD echo
$TRY_CMD echo "Accesses by Date"
$TRY_CMD echo "Accesses by Date" | sed 's/./-/g'
$TRY_CMD uniq -c < "$file_dates"

# Accesses by day of week
$TRY_CMD echo
$TRY_CMD echo "Accesses by Day of Week"
$TRY_CMD echo "Accesses by Day of Week" | sed 's/./-/g'
$TRY_CMD sed 's|/|-|g' "$file_dates" | date -f - +%a 2>/dev/null | sort | uniq -c | sort -rn

# Accesses by Local Hour
$TRY_CMD echo
$TRY_CMD echo "Accesses by Local Hour"
$TRY_CMD echo "Accesses by Local Hour" | sed 's/./-/g'
$TRY_CMD awk -F: '{print $2}' "$file_times" | sort | uniq -c
