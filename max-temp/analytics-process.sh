#!/bin/bash

## Processing
cat "${data_file}" |
  cut -c 88-92 |
  grep -v 999 |
  sort -rn |
  head -n1 > ${outputs_dir}/max.txt

cat "${data_file}" |
  cut -c 88-92 |
  grep -v 999 |
  sort -n |
  head -n1 > ${outputs_dir}/min.txt

cat "${data_file}" |
  cut -c 88-92 |
  grep -v 999 |
  awk "{ total += \$1; count++ } END { print total/count }" > ${outputs_dir}/average.txt 