#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/ml"
input_dir="${eval_dir}/inputs"
parsed_args=()
size="full"
for arg in "$@"; do
    case "$arg" in
        --small)
            parsed_args+=("$arg")
            size="small"
            ;;
        --min)
            parsed_args+=("$arg")
            size="min"
            ;;
    esac
done

mkdir -p $eval_dir/outputs
mkdir -p $eval_dir/inputs/input_"$size"

export TMP="$eval_dir/inputs/input_$size"
# Generating model & samples
python3 $eval_dir/scripts/gen_model.py 100
python3 $eval_dir/scripts/gen_samples.py "${parsed_args[@]}"

URL='https://atlas.cs.brown.edu/data'

full_dir="${input_dir}/images_full"
small_dir="${input_dir}/images_small"
min_dir="${input_dir}/images_min"
models_dir="${input_dir}/models"

if [ ! -d "$models_dir" ]; then
  mkdir -p "$models_dir"
  wget --no-check-certificate "${URL}/models.zip" -O "${input_dir}/models.zip"
  unzip -q "${input_dir}/models.zip" -d "${input_dir}/tmp_models"
  mv "${input_dir}/tmp_models"/models/* "$models_dir"
  rm -r "${input_dir}/tmp_models" "${input_dir}/models.zip"
fi

if [[ "$size" == "min" ]]; then
  if [ -d "$min_dir" ]; then
    echo "Data already downloaded and extracted."
    exit 0
  fi
  mkdir -p "$min_dir"
  cp "${eval_dir}"/min_inputs/* "$min_dir/"
  exit 0
fi
# if size == small
if [[ "$size" == "small" ]]; then
  if [ -d "$small_dir" ]; then
    echo "Data already downloaded and extracted."
    exit 0
  fi
  mkdir -p "$small_dir"
  wget --no-check-certificate "${URL}/pl-06-P_F-A_N-20250401T083751Z-001.zip" -O "${input_dir}/small.zip"
  unzip -q "${input_dir}/small.zip" -d "${input_dir}/tmp_small"
  mv "${input_dir}/tmp_small"/*/* "$small_dir"
  rm -r "${input_dir}/tmp_small" "${input_dir}/small.zip"
  exit 0
fi

if [ -d "$full_dir" ]; then
  echo "Data already downloaded and extracted."
  exit 0
fi
mkdir -p "$full_dir"
wget --no-check-certificate "${URL}/pl-01-PFW-20250401T083800Z-001.zip" -O "${input_dir}/full.zip"
unzip -q "${input_dir}/full.zip" -d "${input_dir}/tmp_full"
mv "${input_dir}/tmp_full"/*/* "$full_dir"
rm -r "${input_dir}/tmp_full" "${input_dir}/full.zip"
