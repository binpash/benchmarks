#!/bin/bash

inputs_dir="$(dirname "$0")/inputs"

mkdir -p "$inputs_dir"

fallocate -l 1G "$inputs_dir/1G.zeros" 
