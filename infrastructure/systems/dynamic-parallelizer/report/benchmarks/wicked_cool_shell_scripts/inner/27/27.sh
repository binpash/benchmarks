#!/bin/sh
IFS='
'

# numberlines--A simple alternative to cat -n, etc.
for filename in "$@"; do
    linecount="1"
    for line in $(cat "$filename"); do
        echo "${linecount}: $line"
        linecount="$((linecount + 1))"
    done
done
