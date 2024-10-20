#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/unix50"
input_dir="${eval_dir}/input"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
scripts_dir="${eval_dir}/scripts"

txt_inputs=$input_dir/full
suffix=.full
size_suffix=_3G.txt
if [[ "$@" == *"--small"* ]]; then
    txt_inputs=$input_dir/small
    suffix=.small
    size_suffix=_1M.txt
fi
txt_outputs=$results_dir/result$suffix
mkdir -p $txt_outputs

inputs=(
    "1.sh;1"
    "2.sh;1"
    "3.sh;1"
    "4.sh;1"
    "5.sh;2"
    "6.sh;3"
    "7.sh;4"
    "8.sh;4"
    "9.sh;4"
    "10.sh;4"
    "11.sh;4"
    "12.sh;4"
    "13.sh;5"
    "14.sh;6"
    "15.sh;7"
    "16.sh;7"
    "17.sh;7"
    "18.sh;8"
    "19.sh;8"
    "20.sh;8"
    "21.sh;8"
    # "22.sh;8"
    "23.sh;9.1"
    "24.sh;9.2"
    "25.sh;9.3"
    "26.sh;9.4"
    # "27.sh;9.5"
    "28.sh;9.6"
    "29.sh;9.7"
    "30.sh;9.8"
    "31.sh;9.9"
    "32.sh;10"
    "33.sh;10"
    "34.sh;10"
    "35.sh;11"
    "36.sh;11"
)
for input in "${inputs[@]}"; do
    IFS=';' read -r script input <<< "$input"
    $scripts_dir/$script $txt_inputs/$input$size_suffix > $txt_outputs/$script.out
done
