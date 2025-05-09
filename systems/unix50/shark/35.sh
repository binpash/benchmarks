#!/bin/bash

# 11.1: year Ritchie and Thompson receive the Hamming medal
# cat $1 | grep 'UNIX' | cut -f 1

grep 'UNIX' < $1 | cut -f 1