#!/bin/bash

[[ -n "$input_file" ]] || echo "script was not provided with \$input_file"
[[ -n "$results_dir" ]] || echo "script was not provided with \$results_dir"

cat "${input_file}" |
  cut -c 88-92 |
  grep -v 999 |
  sort -rn |
  head -n1 > ${results_dir}/max.txt

cat "${input_file}" |
  cut -c 88-92 |
  grep -v 999 |
  sort -n |
  head -n1 > ${results_dir}/min.txt

cat "${input_file}" |
  cut -c 88-92 |
  grep -v 999 |
  awk "{ total += \$1; count++ } END { print total/count }" > ${results_dir}/average.txt 
