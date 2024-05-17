#!/bin/bash
cat ${IN}apache.log  | grep -o "from [^ ]*" | cut -d ' ' -f2 | sort | uniq -c | sort -nr
