#!/bin/bash

# 11.2: most repeated first name in the list?
# cat $1 | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d

cut -f 2 < $1 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d