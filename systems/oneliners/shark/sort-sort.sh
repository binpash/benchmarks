#!/bin/bash
# Calculate sort twice

# cat $1 | tr A-Z a-z | sort | sort -r

tr A-Z a-z < $1 | sort | sort -r