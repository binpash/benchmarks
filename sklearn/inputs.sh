#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p tmp
mkdir -p result

# Currently just dumped the entire dataset, but ideally we actually download it