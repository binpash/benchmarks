#!/bin/bash

sudo apt update 

if ! dpkg -s samtools &> /dev/null; then
    sudo apt install samtools -y
else
    exit 0
fi

