#!/bin/bash

# 8.4: find longest words without hyphens
# cat $1 | tr -c "[a-z][A-Z]" '\n' | sort | awk "length >= 16"

tr -c "[a-z][A-Z]" '\n' < $1 | sort | awk "length >= 16"