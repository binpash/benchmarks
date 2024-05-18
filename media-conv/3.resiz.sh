#!/bin/bash

# tag: resize image 
find $IN -name "*.jpg" | 
  xargs -n1 basename |
  sed "s;\(.*\);-resize 70% $IN/\1 $OUT/\1.70;" |
  xargs -L1  convert
