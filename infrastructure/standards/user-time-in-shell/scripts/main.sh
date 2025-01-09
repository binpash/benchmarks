#!/bin/bash

start=$(date +%s%3N)
for ((i=0; i < 1000000; i++)) {
    :
}
end=$(date +%s%3N)
# want 10% of the time to be in the shell, so sleep the remaining time
duration=$(awk "BEGIN {print 9 * ($end - $start) / 1000}")
sleep $duration

