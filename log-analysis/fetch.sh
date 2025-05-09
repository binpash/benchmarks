#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/inputs"
URL=https://atlas.cs.brown.edu/data
mkdir -p "$input_dir"

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done


if [[ ! -d "$input_dir/nginx-logs_$size" ]] || [[ ! -d "$input_dir/pcaps_$size" ]]; then
    if [[ "$size" == "min" ]]; then
        cp -r "${eval_dir}/min_inputs/nginx-logs" "$input_dir/nginx-logs_$size"
        cp -r "${eval_dir}/min_inputs/pcaps" "$input_dir/pcaps_$size"
        exit 0
    fi

    if [[ "$size" == "small" ]]; then
        wget --no-check-certificate $URL/pcaps.zip -O "$input_dir/pcaps_$size.zip"
        unzip "$input_dir/pcaps_$size.zip" -d "$input_dir"
        mv "$input_dir/pcaps" "$input_dir/pcaps_$size"
        rm "$input_dir/pcaps_$size.zip"

        zip_dst="$input_dir/nginx.zip"
        wget --no-check-certificate $URL/nginx.zip -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        mv "$input_dir/nginx-logs" "$input_dir/nginx-logs_$size"
        rm "$zip_dst"
        exit 0
    fi

    zip_dst="$input_dir/pcaps.zip"
    wget --no-check-certificate "$URL"/pcaps.zip -O "$zip_dst"
    unzip "$zip_dst" -d "$input_dir"
    mv "$input_dir/pcaps/" "$input_dir/pcaps_$size/"
    rm "$zip_dst"

    wget --no-check-certificate "$URL"/pcaps_large.zip -O "$zip_dst"
    unzip "$zip_dst" -d "$input_dir"/pcaps_$size
    rm "$zip_dst"

    zip_dst="$input_dir/nginx.zip"
    wget --no-check-certificate $URL/nginx.zip -O "$zip_dst"
    unzip "$zip_dst" -d "$input_dir"
    mv "$input_dir/nginx-logs" "$input_dir/nginx-logs_$size"
    rm "$zip_dst"

    zip_dst="$input_dir/nginx_large.zip"
    wget --no-check-certificate "$URL"/log-analysis/web-server-access-logs.zip -O "$zip_dst"
    unzip "$zip_dst" -d "$input_dir"
    mv "$input_dir/access.log" "$input_dir/nginx-logs_$size/access.log"
    rm "$zip_dst"
else
    echo "Data already downloaded and extracted."
    exit 0
fi
