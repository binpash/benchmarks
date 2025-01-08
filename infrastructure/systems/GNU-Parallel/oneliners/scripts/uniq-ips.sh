#!/bin/bash

# Finding unique ips 
# from: https://blog.cloudflare.com/when-bloom-filters-dont-bloom/
cat "$1" | sort | uniq

# No opportunity for GNU parallel