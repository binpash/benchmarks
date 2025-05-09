#!/bin/bash

# 5.1: extract hello world
# cat $1 | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12

grep 'print' < $1 | cut -d "\"" -f 2 | cut -c 1-12