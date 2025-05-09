
#!/bin/bash

## Initialize the necessary temporary files
REPO_DIR="$REPO_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
file1="$OUTPUT_DIR/file1.txt"
file2="$OUTPUT_DIR/file2.txt"

# Create list of files
find "$REPO_DIR" -type f |
xargs openssl md5 |
sed 's/^MD5(//;s/)= / /' |
sort -k2 > "$file1"
awk '{print $2}' < "$file1" | uniq -d > "$file2"
join -2 2 "$file2" "$file1" |
awk '
BEGIN {ORS=""}
$1 != prev && prev {print "\n"}
END {if (prev) print "\n"}
{if (prev) print " "; prev = $1; print $2}'
