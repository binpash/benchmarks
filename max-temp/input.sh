#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

URL='https://www1.ncdc.noaa.gov/pub/data/noaa/'
FROM=${FROM:-2015}
TO=${TO:-2015}
sample_starting_index=1234
sample_count=250

mkdir -p "$input_dir"

seq $FROM $TO |
  sed "s;^;$URL;" |
  sed 's;$;/;' |
  xargs -n1 -r curl --insecure |
  grep gz |
  sed "s;.*\"\(.*\)\(20[0-9][0-9]\).gz\".*;$URL\2/\1\2.gz;" |
  tail -n +$sample_starting_index |
  head -n $sample_count |
  xargs -n1 curl --insecure |
  gunzip > "$input_dir/temperatures.full.txt"

head -n 200 "$input_dir/temperatures.full.txt" \
    > "$input_dir/temperatures.small.txt"
