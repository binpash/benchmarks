#!/bin/bash

## This is used to test against the optimal parallel ordering (theoretical best)

for year in $(seq $FROM $TO); do
    (
        echo "Processing year: $year"
        find "$RESOURCE_DIR/$year" -type f -name '*.gz' -print0 | xargs -0 -I{} sh -c "gunzip -c "{}" >> \"$RESOURCE_DIR/$year.txt\""

        # Process in parallel
        cut -c 89-92 "$RESOURCE_DIR/$year.txt" | grep -v 999 | sort -rn | head -n1 &
        cut -c 89-92 "$RESOURCE_DIR/$year.txt" | grep -v 999 | sort -n | head -n1 &
        cut -c 89-92 "$RESOURCE_DIR/$year.txt" | grep -v 999 | awk '{ total += $1; count++ } END { print total/count }' &
        wait
    ) &
done
wait
