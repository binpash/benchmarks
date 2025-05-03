#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
input_dir="${eval_dir}/inputs"

URL='https://atlas.cs.brown.edu/data'
URL=$URL/max-temp/noaa/
FROM=2000
TO=2024

n_samples=99999
suffix="full"

mkdir -p "${input_dir}"

for arg in "$@"; do
    if [[ "$arg" == "--min" ]]; then
      min_inputs="$eval_dir/min_inputs/"
      mkdir -p "$input_dir"
      cp -r "$min_inputs"/* "$input_dir/"
      exit 0
    fi
    if [[ "$arg" == "--small" ]]; then
      FROM=2000
      TO=2000
      n_samples=700
      suffix="small"
    fi
done

seq "$FROM" "$TO" |
  sed "s;^;$URL;" |
  sed 's;$;/;' |
  xargs -n1 -r curl --insecure |
  grep gz |
  sed "s;.*\"\(.*\)\(20[0-9][0-9]\).gz\".*;$URL\2/\1\2.gz;" |
  head -n "$n_samples" |
  xargs -n1 curl --insecure |
  gunzip >"$input_dir/temperatures.$suffix.txt"
