#!/bin/bash
# 8.2: find Bell Labs location where Dennis Ritchie had his office
cat $IN | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2 | awk "{\$1=\$1};1"
