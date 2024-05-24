#!/bin/bash

# 10.2: list Turing award recipients while working at Bell Labs
cat $1 | sed 1d | grep 'Bell' | cut -f 2
