#!/bin/bash

# 4.2: find pieces captured by Belle
cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l
