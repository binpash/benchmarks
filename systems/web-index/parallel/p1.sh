#!/bin/bash
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
WIKI=${WIKI:-$PASH_TOP/web-index}

export WIKI
# Squash all HTML for each URL into a single line, streaming fashion
# It also prefixes with the URL

page_per_line () {
  cat "$WIKI/$0" | tr -d "\n\r" | tr -d '\n' | sed -e '/.$/a\'
}

export -f page_per_line

# xargs:
# add `-t` for debugging
cat $WIKI/input/index.txt | xargs -0 -d '\n' -n 1 bash -c 'page_per_line "$@"'

