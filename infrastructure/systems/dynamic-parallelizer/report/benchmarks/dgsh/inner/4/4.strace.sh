
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

## Initialize the necessary temporary files
REPO_DIR="$REPO_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
file1="$OUTPUT_DIR/file1.txt"
file2="$OUTPUT_DIR/file2.txt"

# Create list of files
logfile=$(generate_unique_file)
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i find "$REPO_DIR" -type f |
xargs openssl md5 |
sed 's/^MD5(//;s/)= / /' |
sort -k2 > "$file1"
awk '{print $2}' < "$file1" | uniq -d > "$file2"


logfile=$(generate_unique_file)
cmd='
BEGIN {ORS=""}
$1 != prev && prev {print "\n"}
END {if (prev) print "\n"}
{if (prev) print " "; prev = $1; print $2}'
strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i join -2 2 "$file2" "$file1" | awk "$cmd" 
