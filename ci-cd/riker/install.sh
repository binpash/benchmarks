#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
scripts_dir="${REPO_TOP}/riker/scripts"

tz="America/New_York"
echo "$tz" | sudo tee /etc/timezone > /dev/null
sudo rm -f /etc/localtime
sudo ln -s "/usr/share/zoneinfo/$tz" /etc/localtime

sudo apt-get update
sudo apt-get install -y build-essential

small_benchmark=(
    "lua"
    "memcached"
    "redis"
    "sqlite"
    "vim"
    "xz"
    "xz-clang"
)

run_small=false

for arg in "$@"; do
    if [ "$arg" = "--small" ]; then
        run_small=true
        break
    fi
done

if [ "$run_small" = true ]; then
    for bench in "${small_benchmark[@]}"; do
        script_path="$scripts_dir/$bench/install.sh"
        if [ -x "$script_path" ]; then
            "$script_path" "$@"
        else
            echo "Error: $script_path not found or not executable."
            exit 1
        fi
    done
    exit 0
fi

for bench in "$scripts_dir"/*; do
    "$bench/install.sh" "$@"
done
