#!/bin/bash

# 4.1: find number of rounds
cat $1 | tr ' ' '\n' | grep '\.' | wc -l
