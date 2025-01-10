#!/bin/bash

# 7.1: identify number of AT&T unix versions
# cat $1 | cut -f 1 | grep 'AT&T' | wc -l

cut -f 1 < $1 | grep 'AT&T' | wc -l