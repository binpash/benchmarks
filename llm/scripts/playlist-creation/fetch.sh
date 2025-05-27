#!/bin/bash
KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
URL='https://atlas.cs.brown.edu/data'
input_dir="$REPO_TOP/llm/inputs/scripts/playlist-creation/inputs"

mkdir -p "$input_dir"

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=small ;;
    --min) size=min ;;
    esac
done

if [[ "$size" == "small" ]]; then
    if [[ -d "$input_dir/songs.small" ]]; then
        echo "Data already downloaded and extracted."
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
    if [[ -d "$input_dir/songs.min" ]]; then
        echo "Data already downloaded and extracted."
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
    if [[ -d "$input_dir/songs.full" ]]; then
        echo "Data already downloaded and extracted."
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
