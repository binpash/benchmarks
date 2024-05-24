#!/bin/bash
# Calculate sort twice

cat $1 | tr A-Z a-z | sort | sort -r
