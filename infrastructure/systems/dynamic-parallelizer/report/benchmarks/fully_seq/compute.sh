#!/bin/bash

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

print_fibonacci "$1"

echo "$1" > "$2"
