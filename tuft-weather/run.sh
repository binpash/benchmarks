#!/usr/bin/env bash

input="input"

# Remove missing values, sort by date and combine month and day into one column
cat "$input" \
    | grep -v "\-99" \
    | sort -n \
    | awk '{ printf "%02d-%02d %s %s\n", $1, $2, $3, $4 }' > formatted.txt

input="formatted.txt"

# Get max and min per month date across all years
cat "$input" \
    | awk '{
        key = sprintf("%s", $1);
        if (!(key in max) || $3 > max[key]) max[key] = $3;
        if (!(key in min) || $3 < min[key]) min[key] = $3;
    } 
    END {
        for (key in max) {
            printf "%s %s %s\n", key, max[key], min[key];
        }
    }' \
    | sort -n > max_min.txt

# For each year, get temperature per day and plot
for year in $(seq 1995 2013); do
    join -1 1 -2 1 -o "0 1.2 1.3 2.2" -e "NA" max_min.txt <(grep "$year" "$input" | cut -d' ' -f 1,3) \
    | grep -v "NA" \
    | ./scripts/plot.py "$year" > "$year.png"
done
