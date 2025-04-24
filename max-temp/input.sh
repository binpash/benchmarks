#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
input_dir="${eval_dir}/inputs"

URL='https://www1.ncdc.noaa.gov/pub/data/noaa/'
FROM=${FROM:-2015}
TO=${TO:-2015}
sample_starting_index=1234
sample_count=250

mkdir -p "${input_dir}"

for arg in "$@"; do
    if [[ "$arg" == "--min" ]]; then
      min_inputs="$eval_dir/min_inputs/"
      mkdir -p "$input_dir"
      cp -r "$min_inputs"/* "$input_dir/"
      exit 0
    fi
done

if [[ -d "$input_dir" ]]; then
  echo "Data already downloaded and extracted."
  exit 0
fi

seq "$FROM" "$TO" |
  sed "s;^;$URL;" |
  sed 's;$;/;' |
  xargs -n1 -r curl --insecure |
  grep gz |
  sed "s;.*\"\(.*\)\(20[0-9][0-9]\).gz\".*;$URL\2/\1\2.gz;" |
  tail -n +$sample_starting_index |
  head -n $sample_count |
  xargs -n1 curl --insecure |
  gunzip >"$input_dir/temperatures.full.txt"

head -n 200 "$input_dir/temperatures.full.txt" \
  >"$input_dir/temperatures.small.txt"

head -n 20 "$input_dir/temperatures.full.txt" \
  >"$input_dir/temperatures.min.txt"
