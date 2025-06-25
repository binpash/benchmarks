#!/bin/bash

# 6.1: order the bodies by how easy it would be to land on them in Thompson's Space Travel game when playing at the highest simulation scale
cat $1 | awk "{print \$2, \$0}" | sort -nr | cut -d ' ' -f 2
