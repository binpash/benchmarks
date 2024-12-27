#!/bin/bash

sudo apt update 

sudo apt install -y --no-install-recommends \
    bash \
    curl \
    grep \
    awk \
    iptables \
    ufw \
    systemd
