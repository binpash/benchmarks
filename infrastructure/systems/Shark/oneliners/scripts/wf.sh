#!/bin/bash
# Calculate the frequency of each word in the document, and sort by frequency

# cat $1 | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | tr A-Z a-z | sort | uniq -c | sort -rn 

tr -c 'A-Za-z' '[\n*]' < $1 | grep -v "^\s*$" | tr A-Z a-z | sort | uniq -c | sort -rn