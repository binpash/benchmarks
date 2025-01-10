#!/bin/bash

# 1.2: extract names and sort
# cat $1 | head -n 2 | cut -d ' ' -f 2

head -n 2 < $1 | cut -d ' ' -f 2