#!/bin/bash

# 8.1: count unix birth-year
# cat $1 | tr ' ' '\n' | grep 1969 | wc -l

tr ' ' '\n' < $1 | grep 1969 | wc -l