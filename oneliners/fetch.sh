#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
URL="https://atlas.cs.brown.edu/data"

input_dir="${TOP}/oneliners/inputs"
mkdir -p "$input_dir"
cd "$input_dir" || exit 1

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done


if [ ! -f ./1M.txt ]; then
    wget --no-check-certificate "$URL"/dummy/1M.txt
    # Add newline to file
    echo >>1M.txt
    dos2unix 1M.txt
fi

if [ ! -f ./dict.txt ]; then
    wget -O - "$URL"/dummy/dict.txt --no-check-certificate | sort >dict.txt
fi

if [ ! -f ./all_cmds.txt ]; then
    wget -O - "$URL"/dummy/all_cmds.txt --no-check-certificate >all_cmds.txt
fi

if [ ! -f ./all_cmdsx100.txt ]; then
    touch all_cmdsx100.txt
    for ((i = 0; i < 100; i++)); do
        cat all_cmds.txt >>all_cmdsx100.txt
    done
fi

# For uniq-ips
if [ "$size" = "small" ]; then
    N=400000 # 4K
elif [ "$size" = "min" ]; then
    N=40
else
    N=40000000 # 40M
fi

../scripts/gen_ips.py "$N" >logs-popcount-org.txt

if [[ "$size" == "small" ]]; then
    if [ ! -f ./10M.txt ]; then
        touch 10M.txt
        for ((i = 0; i < 10; i++)); do
            cat 1M.txt >>10M.txt
        done
    fi

    if [ ! -f ./30M.txt ]; then
        touch 30M.txt
        for ((i = 0; i < 3; i++)); do
            cat 10M.txt >>30M.txt
        done
    fi
elif [[ "$size" == "min" ]]; then
    exit 0
fi

# full size files
if [ ! -f ./1G.txt ]; then
    touch 1G.txt
    for ((i = 0; i < 1000; i++)); do
        cat 1M.txt >>1G.txt
    done
fi

if [ ! -f ./3G.txt ]; then
    touch 3G.txt
    for ((i = 0; i < 3; i++)); do
        cat 1G.txt >>3G.txt
    done
fi
