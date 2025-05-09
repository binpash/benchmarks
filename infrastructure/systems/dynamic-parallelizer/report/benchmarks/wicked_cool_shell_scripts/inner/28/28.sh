#!/bin/sh
# toolong--Feeds the fmt command only those lines in the input stream
# that are longer than the specified length
IFS='
'

width=72
if [ ! -r "$1" ] ; then
    echo "Cannot read file $1" >&2
    echo "Usage: $0 filename" >&2
    exit 1
fi
for input in $(cat "$1")
do
    varlength="$(echo "$input" | wc -c | sed 's/[^[:digit:]]//g')"
    if [ $varlength -gt $width ] ; then
        echo "$input" | fmt
    else
        echo "$input"
    fi
done
