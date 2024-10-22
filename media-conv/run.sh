#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
input_dir="${eval_dir}/input"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
mkdir -p $results_dir

img_convert_input="$input_dir/jpg_full/jpg"
to_mp3_input="$input_dir/wav_full"
suffix=".full"

if [[ "$@" == *"--small"* ]]; then
    img_convert_input="$input_dir/jpg_small/jpg"
    to_mp3_input="$input_dir/wav_small"
    suffix=".small"
fi

$scripts_dir/img_convert.sh $img_convert_input $results_dir/img_convert$suffix > $results_dir/img_convert$suffix.log
$scripts_dir/to_mp3.sh $to_mp3_input $results_dir/to_mp3$suffix > $results_dir/to_mp3$suffix.log
