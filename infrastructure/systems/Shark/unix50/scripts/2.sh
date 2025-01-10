#!/bin/bash

# 1.1: extract names and sort
# cat $1 | cut -d ' ' -f 2 | sort

cut -d ' ' -f 2 < $1 | sort