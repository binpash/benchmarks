#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

FROM=${FROM:-2015}
TO=${TO:-2015}
URL='https://www1.ncdc.noaa.gov/pub/data/noaa/'

## Downloading and extracting
seq $FROM $TO |
  sed "s;^;$URL;" |
  sed 's;$;/;' |
  xargs -r -n1 --insecure |
  grep gz |
  tr -s ' \n' |
  cut -d ' ' -f9 |
  sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
  sed "s;^;$URL;" |
  xargs -n1 curl --insecure |
  gunzip > "$input_dir/temperatures2015.txt"
