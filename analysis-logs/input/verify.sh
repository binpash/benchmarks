#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-logs/"
results_dir="${eval_dir}/results"
input_dir="${eval_dir}/input"

if [ "$(md5sum $results_dir/* | awk '{print $1}')" == "$(cat $input_dir/checksum.md5 | awk '{print $1}')" ];
then
    echo "Valid"
else
    echo "Invalid"
    exit 1
fi
