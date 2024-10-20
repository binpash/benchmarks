#!/bin/bash
# compress all files in a directory
mkdir -p $2

for item in $1/*.pcapng;
do
    output_name="$2/$(basename $item).zip"
    cat $item | gzip --no-name -c > $output_name
done
