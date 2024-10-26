#!/bin/bash

# xargs:
# add `-t` for debugging

cat "$INPUT_FILE" | xargs -d '\n' -I {} bash "$SCRIPT_DIR/page_per_line.sh" {} $@ > out