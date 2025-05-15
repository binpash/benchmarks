#!/bin/bash
KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
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
    data_url=https://os.unil.cloud.switch.ch/fma/fma_small.zip
    zip_dst=$input_dir/fma_small.zip
    out_dir=$input_dir
    wget --no-check-certificate $data_url -O $zip_dst || {
        echo "Failed to download $data_url"
        exit 1
    }
    unzip -qq $zip_dst -d $out_dir || {
        echo "Failed to unzip $zip_dst"
        exit 1
    }
    rm "$zip_dst"
    mv "$out_dir" "$input_dir/songs.small/"
    # take out problematic files
    cd "$input_dir/songs.small/" || exit 1
    rm -rf 11 21 29 54 98 99 108 133
    exit 0
elif [[ "$size" == "min" ]]; then
    if [[ -d "$input_dir/songs.min" ]]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi
    cp -r "$REPO_TOP/llm/scripts/playlist-creation/min_inputs/"* "$input_dir/songs.min/"
    exit 0
else
    if [[ -d "$input_dir/songs.full" ]]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi

    echo "Downloading full dataset."
    data_url=https://os.unil.cloud.switch.ch/fma/fma_medium.zip
    zip_dst="$input_dir/fma_medium.zip"
    out_dir="$input_dir"
    wget --no-check-certificate $data_url -O $zip_dst
    unzip -qq $zip_dst -d $out_dir
    rm "$zip_dst"
    mv "$out_dir" "$input_dir/songs.full/"
fi