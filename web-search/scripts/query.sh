#!/bin/bash

# Provided an appropriate index, the query could be implemented using grep
# along  with appropriate processing and stemming of the input strings.

grep "$(echo "$@" | ./c/process.sh | ./c/stem.js | tr "\r\n" "  ")" d/global-index.txt
