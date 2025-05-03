#!/bin/bash

sudo apt update 

sudo apt install -y --no-install-recommends \
    sudo \
    curl \
    wget \
    unzip \
    gzip \
    coreutils \
    ffmpeg unrtf imagemagick zstd \
    git \
    xz-utils
