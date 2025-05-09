#!/bin/bash
# 10.3: extract Ritchie's username
cat $IN | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'
