#!/bin/bash

# 1.0: extract the last name
cat $1 | cut -d ' ' -f 2
