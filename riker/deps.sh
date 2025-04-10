#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"

tz="America/New_York"
echo "$tz" | sudo tee /etc/timezone > /dev/null
sudo rm -f /etc/localtime
sudo ln -s "/usr/share/zoneinfo/$tz" /etc/localtime

sudo apt update

sudo apt install -y build-essential

for bench in "$scripts_dir"/*; do
    "$bench/deps.sh" "$@"
done

