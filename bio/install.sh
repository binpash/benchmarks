#!/bin/bash

sudo apt-get update 

if ! dpkg -s samtools &> /dev/null; then
    sudo apt install samtools -y
fi
