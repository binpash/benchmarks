#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

input_files=("1M.txt" "1G.txt" "3G.txt" "all_cmds.txt" "all_cmdsx100.txt" "dict.txt")

if [ ! -f ./1M.txt ]; then
    wget --no-check-certificate https://atlas-group.cs.brown.edu/data/dummy/1M.txt
    # TODO: Add newline to the original file
    echo >> 1M.txt
    dos2unix 1M.txt
fi

if [ ! -f ./1G.txt ]; then
    touch 1G.txt
    for (( i = 0; i < 1000; i++ )); do
        cat 1M.txt >> 1G.txt
    done
fi

if [ ! -f ./3G.txt ]; then
    touch 3G.txt
    for (( i = 0; i < 3; i++ )); do
        cat 1G.txt >> 3G.txt
    done
fi

if [ ! -f ./dict.txt ]; then
    wget -O - https://atlas-group.cs.brown.edu/data/dummy/dict.txt --no-check-certificate | LC_ALL=C sort > dict.txt
fi

if [ ! -f ./all_cmds.txt ]; then
    # TODO: Upload this file to the server
    cp ../all_cmds.txt all_cmds.txt
fi

if [ ! -f ./all_cmdsx100.txt ]; then
    touch all_cmdsx100.txt
    for (( i = 0; i < 100; i++ )); do
        cat all_cmds.txt >> all_cmdsx100.txt
    done
fi

# For uniq-ips
if [ "$1" = "--small" ]; then
    N=4000 # 4K
else
    N=40000000 # 40M
fi

../scripts/gen_ips.py "$N" > logs-popcount-org.txt
