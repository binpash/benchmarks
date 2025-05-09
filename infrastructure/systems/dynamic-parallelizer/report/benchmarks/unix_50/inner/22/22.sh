#!/bin/bash
# # 8.5: Find second-most-freq 8-character word(s) without hyphens
cat $IN | tr -cs '[:alpha:]' '\n' | awk 'length($0) == 8' | tr '[:upper:]' '[:lower:]' | sort | uniq -c | sort -nr | awk 'NR==2 {print $2}'
