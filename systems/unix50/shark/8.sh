#!/bin/bash

# 4.2: find pieces captured by Belle
# cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

tr ' ' '\n' < $1 | grep 'x' | grep '\.' | wc -l