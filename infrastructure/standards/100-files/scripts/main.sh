#!/bin/bash

run() {
    for f in "$(dirname "$0")/../inputs"/*; do
        # keep the file open for 0.1 seconds
        exec {fd}<> "$f"
        sleep 0.1 
        exec {fd}>&-
    done
}

run
