#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="$REPO_TOP/llm"
input_dir="$eval_dir/inputs"

mkdir -p "$input_dir"

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done
export LC_ALL=C
URL="https://atlas-group.cs.brown.edu/data"

if [[ "$size" == "small" ]]; then
    # if inputs exist
    if [[ -d "$input_dir/jpg.small" ]]; then
        echo "Image data already downloaded and extracted."
    else
        data_url="${URL}"/small/jpg.zip
        zip_dst=$input_dir/jpg.small.zip
        out_dir=$input_dir/jpg.small
        wget --no-check-certificate $data_url -O $zip_dst || {
            echo "Failed to download $data_url"
            exit 1
        }
        unzip $zip_dst -d $out_dir || {
            echo "Failed to unzip $zip_dst"
            exit 1
        }
        rm "$zip_dst"
    fi
        if [[ -d "$input_dir/songs.small" ]]; then
        echo "Song already downloaded and extracted."
        exit 0
    fi
    data_url="${URL}/llm/playlist_small.tar.gz"
    wget --no-check-certificate $data_url -O "$input_dir"/playlist_small.tar.gz || {
        echo "Failed to download $data_url"
        exit 1
    }
    tar -xzf "$input_dir/playlist_small.tar.gz" -C "$input_dir" || {
        echo "Failed to extract $input_dir/playlist_small.tar.gz"
        exit 1
    }
    rm "$input_dir/playlist_small.tar.gz"
    mv "$input_dir/playlist_small" "$input_dir/songs.small"
    exit 0

elif [[ "$size" == "min" ]]; then
    if [[ -d "$input_dir/jpg.min" ]]; then
        echo "Image data already downloaded and extracted."
    else
        min_inputs="$eval_dir/min_inputs/"
        out_dir="$input_dir/jpg.min/jpg"
        mkdir -p "$out_dir"
        cp -r "$min_inputs"/* "$out_dir/"
    fi
        if [[ -d "$input_dir/songs.min" ]]; then
        echo "Song data already downloaded and extracted."
        exit 0
    fi
    data_url="${URL}/llm/playlist_min.tar.gz"
    wget --no-check-certificate $data_url -O $input_dir/playlist_min.tar.gz || {
        echo "Failed to download $data_url"
        exit 1
    }
    tar -xzf "$input_dir/playlist_min.tar.gz" -C "$input_dir" || {
        echo "Failed to extract $input_dir/playlist_min.tar.gz"
        exit 1
    }
    rm "$input_dir/playlist_min.tar.gz"
    mv "$input_dir/playlist_min" "$input_dir/songs.min"
    exit 0
else
    if [[ -d "$input_dir/jpg" ]]; then
        echo "Image data already downloaded and extracted."
    else
        echo "Downloading full dataset."
        data_url=https://atlas-group.cs.brown.edu/data/full/jpg.zip
        zip_dst="$input_dir/jpg.zip"
        out_dir="$input_dir/jpg"
        wget --no-check-certificate $data_url -O $zip_dst
        unzip $zip_dst -d $out_dir
        rm "$zip_dst"
    fi
        if [[ -d "$input_dir/songs.full" ]]; then
        echo "Song data already downloaded and extracted."
        exit 0
    fi
    echo "Downloading full dataset."
    data_url="${URL}/llm/playlist_full.tar.gz"
    wget --no-check-certificate $data_url -O "$input_dir"/playlist_full.tar.gz || {
        echo "Failed to download $data_url"
        exit 1
    }
    tar -xzf "$input_dir/playlist_full.tar.gz" -C "$input_dir" || {
        echo "Failed to extract $input_dir/playlist_full.tar.gz"
        exit 1
    }
    rm "$input_dir/playlist_full.tar.gz"
    mv "$input_dir/playlist_full" "$input_dir/songs.full"
fi