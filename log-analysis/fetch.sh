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
    --min) size=min ;;
    esac
done

if [ -f "$input_dir/routeviews.mrt" ]; then
    echo "routeviews.mrt already exists, skipping download."
else
    wget --no-check-certificate https://archive.routeviews.org/route-views.linx/bgpdata/2019.10/RIBS/rib.20191001.0000.bz2 -O "$input_dir/rib.20191001.0000.bz2"
    bunzip2 "$input_dir/rib.20191001.0000.bz2"
    rm "$input_dir/rib.20191001.0000.bz2"
    mv "$input_dir/rib.20191001.0000" "$input_dir/routeviews.mrt"
fi

if [[ ! -d "$input_dir/nginx-logs_$size" ]] || [[ ! -d "$input_dir/pcaps_$size" ]] || [[ ! -d "$input_dir/port_scan_$size" ]]; then
    if [[ "$size" == "min" ]]; then
        mkdir -p "$input_dir/nginx-logs_$size"
        mkdir -p "$input_dir/pcaps_$size"
        cp "${eval_dir}/min_inputs/nginx-logs/"* "$input_dir/nginx-logs_$size"
        cp "${eval_dir}/min_inputs/pcaps/"* "$input_dir/pcaps_$size"

        mkdir -p "$input_dir/port_scan_$size"
        for log in "$input_dir/nginx-logs_$size/"*; do
            python3 "${eval_dir}/scripts/format_inputs.py" "$log" >"$input_dir/port_scan_$size/$(basename "$log")"
        done
        # concat logs in port_scan_$size
        cat "$input_dir/port_scan_$size/"* >"$input_dir/port_scan_$size/all_logs.jsonl"
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

        mkdir -p "$input_dir/port_scan_$size"
        for log in "$input_dir/nginx-logs_$size"/*; do
            python3 "${eval_dir}/scripts/format_inputs.py" "$log" >"$input_dir/port_scan_$size/$(basename "$log")"
        done
        # concat logs in port_scan_$size
        cat "$input_dir/port_scan_$size/"* >"$input_dir/port_scan_$size/all_logs.jsonl"
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
    mkdir -p "$input_dir/port_scan_$size"
    for log in "$input_dir/nginx-logs_$size"/*; do
        python3 "${eval_dir}/scripts/format_inputs.py" "$log" >"$input_dir/port_scan_$size/$(basename "$log")"
    done
    cat "$input_dir/port_scan_$size/"* >"$input_dir/port_scan_$size/all_logs.jsonl"
else
    echo "Data already downloaded and extracted."
    exit 0
fi
