#!/bin/bash

sudo apt-get update 

pkgs="bash curl grep gawk iptables ufw procps net-tools fail2ban iproute2 git patch time"

for pkg in $pkgs; do
    if ! dpkg -l | grep -q "^ii  $pkg "; then
        sudo apt-get install -y "$pkg"
    fi
done
