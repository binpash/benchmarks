#!/bin/bash

start=$(date +%s%3N)
i=0
while [ "$i" -lt 1000000 ]; do
    i=$((i + 1))
done
end=$(date +%s%3N)
# want 10% of the time to be in the shell, so sleep the remaining time
duration=$(awk "BEGIN {print 9 * ($end - $start) / 1000}")
sleep $duration

