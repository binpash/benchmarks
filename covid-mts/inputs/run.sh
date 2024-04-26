#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/covid-mts"
inputs_dir="${eval_dir}/inputs"
results_dir="${eval_dir}/results"

shell="/bin/bash"

mkdir -p $results_dir

export IN="${inputs_dir}/in.csv"

for i in 1 2 3 4 5
do
    script="${eval_dir}/${i}.sh"
    echo "Executing $script..."
    $shell "$script" > "$results_dir/$i.out"
done
