#!/bin/bash


[[ -n "$input_file" ]] || { echo "script was not provided with \$input_file"; exit 1; }
[[ -n "$statistics_dir" ]] || { echo "script was not provided with \$statistics_dir"; exit 1; }

mkdir -p "${statistics_dir}"

tee < "${input_file}" | cut -c 89-92 | grep -v 999 | sort -rn | head -n1 > "${statistics_dir}/max.txt" &
tee < "${input_file}" | cut -c 89-92 | grep -v 999 | sort -n | head -n1 > "${statistics_dir}/min.txt" &
tee < "${input_file}" | cut -c 89-92 | grep -v 999 | awk '{ total += $1; count++ } END { print total/count }' > "${statistics_dir}/average.txt" &

wait

