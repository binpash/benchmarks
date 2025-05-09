
#!/bin/bash

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)

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
