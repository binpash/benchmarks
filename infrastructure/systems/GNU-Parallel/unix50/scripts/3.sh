#!/bin/bash

# 1.2: extract names and sort
cat $1 | head -n 2 | cut -d ' ' -f 2

# we only edit the first two lines of the file
# No GNU Parallel transformation in this script