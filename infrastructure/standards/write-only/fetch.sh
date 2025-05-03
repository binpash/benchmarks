#!/bin/bash

base="$(dirname "$0")"

head -c 1000000 /dev/zero | tr '\0' '1' > "$base/1M.ones"
