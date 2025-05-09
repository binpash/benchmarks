#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file="$1"

if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 2
fi

print_fibonacci() {
    local n="$1"
    local a=0
    local b=1
    
    for (( i=0; i<n; i++ )); do
        # echo -n "$a "
        local temp="$a"
        a="$b"
        b=$((temp + b))
    done
    # echo
}

input_file_number=$(head -n 1 "$input_file")

print_fibonacci "$input_file_number"

# Trigger a write dependency
exec 3<>"$input_file"
