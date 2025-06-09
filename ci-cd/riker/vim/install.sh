#!/bin/bash

sudo apt-get update
sudo apt-get install -y --no-install-recommends git \
    gcc \
    make \
    libncurses-dev \
    libsm-dev \
    libice-dev \
    libxt-dev \
    libx11-dev \
    libxdmcp-dev \
    libselinux-dev
