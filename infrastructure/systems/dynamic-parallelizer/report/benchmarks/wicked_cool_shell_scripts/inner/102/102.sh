#!/bin/sh
# bulkrename--Renames specified files by replacing text in the filename
# Usage: $0 <find> <replace>

if [ $# -ne 2 ]; then
    echo "Usage: $0 <find> <replace>"
    exit 1
fi

match=$1
replace=$2

for i in *"$match"*; do
    if [ -e "$i" ]; then
        newname=$(echo "$i" | sed "s/$match/$replace/")
        mv "$i" "$newname"
        echo "Renamed file $i to $newname"
    fi
done
