#!/bin/bash

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

# Consistent sorting
export LC_ALL=C

export TZ=$(date +%Z)

# Print initial header only if DGSH_DRAW_EXIT is not set
if [ -z "${DGSH_DRAW_EXIT}" ]
then
    logfile=$(generate_unique_file)
    strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -e "\t\tWWW server statistics\n\t\t=====================\n\nSummary\n======="
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
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cat "$INPUT_FILE" > "$file_initial"

# Number of accesses
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of accesses: '

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i wc -l < "$file_initial"

# Total transferred bytes
cmd='{s += $NF} END {print s}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" "$file_initial" > "$file_bytes"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of Gbytes transferred: '

cmd='{print $1 / 1024 / 1024 / 1024}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" "$file_bytes"

# Process Host names
cmd='{print $1}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" "$file_initial" > "$file_hosts"

# Number of accesses
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of accesses: '

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i wc -l < "$file_hosts"

# Sorted hosts
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i sort "$file_hosts" > "$file_sorted_hosts"

# Unique hosts
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq "$file_sorted_hosts" > "$file_unique_hosts"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of hosts: '

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i wc -l < "$file_unique_hosts"

# Number of TLDs
cmd='$NF !~ /[0-9]/ {print $NF}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F. "$cmd" "$file_unique_hosts" | sort -u | wc -l

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of top level domains: '

# Top 10 hosts
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 10 Hosts"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 10 Hosts" | sed 's/./-/g'

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq -c "$file_sorted_hosts" | sort -rn | head -10

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

# Top 20 TLDs
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 20 Level Domain Accesses"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 20 Level Domain Accesses" | sed 's/./-/g'

cmd='$NF !~ /^[0-9]/ {print $NF}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F. "$cmd" "$file_sorted_hosts" | sort | uniq -c | sort -rn | head -20

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

# Domains
cmd='BEGIN {OFS = "."} $NF !~ /^[0-9]/ {$1 = ""; print}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F. "$cmd" "$file_sorted_hosts" | sort > "$file_domains"

# Number of domains
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of domains: '

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq "$file_domains" | wc -l

# Top 10 domains
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 10 domains"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 10 domains" | sed 's/./-/g'

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq -c "$file_domains" | sort -rn | head -10

# Hosts by volume
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 10 Hosts by Transfer"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 10 Hosts by Transfer" | sed 's/./-/g'

cmd='{bytes[$1] += $NF} END {for (h in bytes) print bytes[h], h}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" "$file_initial" | sort -rn | head -10

# Sorted page name requests
cmd='{print $7}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" "$file_initial" | sort > "$file_requests"

# Top 20 area requests (input is already sorted)
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 20 area requests"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 20 area requests" | sed 's/./-/g'

cmd='{print $2}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F/ "$cmd" "$file_requests" | uniq -c | sort -rn | head -20

# Number of different pages
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Number of different pages: '

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq "$file_requests" | wc -l

# Top 20 requests
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 20 requests"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Top 20 requests" | sed 's/./-/g'

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq -c "$file_requests" | sort -rn | head -20


# Access time: dd/mmm/yyyy:hh:mm:ss
cmd='{print substr($4, 2)}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk "$cmd" "$file_initial" > "$file_times"

# Just dates
cmd='{print $1}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F: "$cmd" "$file_times" > "$file_dates"

# Number of days
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo -n 'Accesses per day: '

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq "$file_dates" | wc -l > "$file_day_count"

logfile=$(generate_unique_file)
cmd='{print NXBYTES / $1 / 1024 / 1024}'
awk -v NXBYTES=$(<"$file_bytes") "$cmd" "$file_day_count"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Accesses by Date"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Accesses by Date" | sed 's/./-/g'

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i uniq -c < "$file_dates"

# Accesses by day of week
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Accesses by Day of Week"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Accesses by Day of Week" | sed 's/./-/g'

cmd='s|/|-|g'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i sed "$cmd" "$file_dates" | date -f - +%a 2>/dev/null | sort | uniq -c | sort -rn

# Accesses by Local Hour
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Accesses by Local Hour"

logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i echo "Accesses by Local Hour" | sed 's/./-/g'

cmd='{print $2}'
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i awk -F: "$cmd" "$file_times" | sort | uniq -c
