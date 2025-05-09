#!/bin/bash


mkdir -p "$2"


pure_func() {
    tempfile=$(mktemp)

    tee "$tempfile" | cut -d "\"" -f3 | cut -d ' ' -f2 | sort | uniq -c | sort -rn
    awk '{print $9}' "$tempfile" | sort | uniq -c | sort -rn
    awk '($9 ~ /404/)' "$tempfile" | awk '{print $7}' | sort | uniq -c | sort -rn
    awk '($9 ~ /502/)' "$tempfile" | awk '{print $7}' | sort | uniq -c | sort -r
    awk -F\" '($2 ~ "/wp-admin/install.php"){print $1}' "$tempfile" | awk '{print $1}' | sort | uniq -c | sort -r
    awk '($9 ~ /404/)' "$tempfile" | awk -F\" '($2 ~ "^GET .*.php")' | awk '{print $7}' | sort | uniq -c | sort -r | head -n 20
    awk -F\" '{print $2}' "$tempfile" | awk '{print $2}' | sort | uniq -c | sort -r
    awk -F\" '($2 ~ "ref"){print $2}' "$tempfile" | awk '{print $2}' | sort | uniq -c | sort -r

    rm -f "$tempfile"
}
export -f pure_func

for log in "$1"/*; do
    logname="$2/$(basename "$log")"
    pure_func < "$log" > "$logname" &
done
wait
