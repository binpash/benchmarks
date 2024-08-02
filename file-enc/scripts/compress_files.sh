#!/bin/bash
# compress all files in a directory
mkdir -p $2

for item in $(ls $1);
do
    output_name=$(basename $item).zip
    cat $item | gzip -c > $2/$output_name
done

echo 'done';
