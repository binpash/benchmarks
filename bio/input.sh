#!/bin/bash

IN=inputs
IN_NAME=input.txt

if [[ "$@" == *"--small"* ]]; then
    IN_NAME=input_small.txt
fi

if [[ $1 == "-c" ]]; then
    rm -rf *.bam
    rm -rf *.sam
    rm -rf ../output
    exit
fi

cd "$(realpath $(dirname "$0"))"

mkdir -p inputs
mkdir -p outputs

cat ${IN_NAME} | while read s_line;
	do
    sample=$(echo $s_line |cut -d " " -f 2);
    if [[ ! -f "inputs/$sample".bam ]]; then
        pop=$(echo $s_line |cut -f 1 -d " ");
        link=$(echo $s_line |cut -f 3 -d " ");
        wget -O "${IN}/$sample".bam "$link"
    fi
done;
