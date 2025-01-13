#!/bin/bash

while read -r line
do
    iconv -c -t ascii//TRANSLIT < "$line" |
    pandoc +RTS -K64m -RTS --from html --to plain --quiet
done
