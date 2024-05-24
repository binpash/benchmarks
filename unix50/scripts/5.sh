#!/bin/bash

# 2.1: get all Unix utilities
cat $1 | cut -d ' ' -f 4 | tr -d ','
