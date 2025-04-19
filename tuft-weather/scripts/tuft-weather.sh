#!/usr/bin/env bash

input="$1"
shift
city="Athens"

# Remove missing values, sort by date and combine month and day into one column
cat "$input" \
    | grep "$city" \
    | grep -v "\-99" \
    | awk '{ printf "%02d-%02d %s %s\n", $1, $2, $3, $4 }' \
    | sort -n > formatted.txt

input="formatted.txt"

# Get max and min per month date across all years
cat "$input" \
    | awk '{
        key = sprintf("%s", $1);
        count[key]++;
        sum[key] += $3;
        sum_sq[key] += $3 * $3;
        if (!(key in max) || $3 > max[key]) max[key] = $3;
        if (!(key in min) || $3 < min[key]) min[key] = $3;
    } 
    END {
        for (key in max) {
            mean = sum[key] / count[key];
            variance = (sum_sq[key] / count[key]) - (mean * mean);
            stddev = (variance > 0) ? sqrt(variance) : 0;
            confidence_delta = 1.96 * stddev / sqrt(count[key]);
            normal_range_low = mean - confidence_delta;
            normal_range_high = mean + confidence_delta;
            printf "%s %s %s %.2f %.2f\n", key, min[key], max[key], normal_range_low, normal_range_high;
        }
    }' \
    | sort -n > processed.txt

mkdir -p "plots/$city"

START_YEAR=$(cat "$input" | head -n 1 | cut -d' ' -f 2)
END_YEAR=$(cat "$input" | tail -n 1 | cut -d' ' -f 2)

# For each year, get temperature per day and plot
for year in $(seq $START_YEAR $END_YEAR); do
    grep "$year" "$input" | cut -d' ' -f 1,3 > "$year.txt"
    comm -12 "$year.txt" <(cut -d' ' -f 1,3 processed.txt) > "max_$year.txt"
    comm -12 "$year.txt" <(cut -d' ' -f 1,2 processed.txt) > "min_$year.txt"

    join -a 1 -a 2 -e "-99" -o "0,1.2,2.2" "$year.txt" "min_$year.txt" > "joined_$year.tmp.txt"
    join -a 1 -a 2 -e "-99" -o "0,1.2,1.3,2.2" "joined_$year.tmp.txt" "max_$year.txt" > "joined_$year.txt"

    join -1 1 -2 1 -e "NA" processed.txt "joined_$year.txt" \
    | grep -v "NA" \
    | sed 's/-99/NaN/g' \
    | ./plot.py "$year" "$city" > "plots/$city/$year.png"
done