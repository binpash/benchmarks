#!/bin/bash

# 9.7: Four corners
# cat $1 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'

sed 2d < $1 | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'