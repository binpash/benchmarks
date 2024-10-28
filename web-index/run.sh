#!/bin/bash

cd "$(dirname "$0")"

./generate_index.sh articles

./web-index/run