#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
input_dir="${eval_dir}/input"
mkdir -p "$input_dir"

small=false
for arg in "$@"; do
  case $arg in
    --small)
      small=true
      shift
      ;;
  esac
done

full_dir="${input_dir}/images_full"
small_dir="${input_dir}/images_small"
models_dir="${input_dir}/models"
mkdir -p "$full_dir" "$small_dir" "$models_dir"

wget --no-check-certificate "https://atlas.cs.brown.edu/data/models.zip" -O "${input_dir}/models.zip"
unzip -q "${input_dir}/models.zip" -d "${input_dir}/tmp_models"
mv "${input_dir}/tmp_models"/models/* "$models_dir"
rm -r "${input_dir}/tmp_models" "${input_dir}/models.zip"

if $small; then
    wget --no-check-certificate "https://atlas.cs.brown.edu/data/pl-06-P_F-A_N-20250401T083751Z-001.zip" -O "${input_dir}/small.zip"
    unzip -q "${input_dir}/small.zip" -d "${input_dir}/tmp_small"
    mv "${input_dir}/tmp_small"/*/* "$small_dir"
    rm -r "${input_dir}/tmp_small" "${input_dir}/small.zip"
    exit 0
fi

wget --no-check-certificate "https://atlas.cs.brown.edu/data/pl-01-PFW-20250401T083800Z-001.zip" -O "${input_dir}/full.zip"
unzip -q "${input_dir}/full.zip" -d "${input_dir}/tmp_full"
mv "${input_dir}/tmp_full"/*/* "$full_dir"
rm -r "${input_dir}/tmp_full" "${input_dir}/full.zip"
