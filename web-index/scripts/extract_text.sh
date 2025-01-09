#!/bin/bash

while read -r line
do
    cat $line |
        iconv -c -t ascii//TRANSLIT |
        pandoc +RTS -K64m -RTS --from html --to plain --quiet
done
