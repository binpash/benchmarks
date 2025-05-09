#!/bin/bash
# # 9.5: backwards running clock, in a backwards poem
cat $IN | rev | cut -c -3 | rev | grep -o '[A-Z][A-Z]*' | rev | tr -d '\n' | rev
