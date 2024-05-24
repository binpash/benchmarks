#!/bin/bash

# 1.3: sort top first names
cat $1 | cut -d ' ' -f 1 | sort | uniq -c | sort -r
