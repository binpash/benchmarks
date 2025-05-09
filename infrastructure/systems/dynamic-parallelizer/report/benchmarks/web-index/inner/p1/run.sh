#!/bin/bash

for line in $(cat "$INPUT_FILE"); do
    cat "$WIKI/$line" | tr -d "\n\r"
done
