#!/bin/bash

cat "$1" | sort | uniq > "out.txt"
