#!/bin/bash

echo "$1" >>${OUT}/visited.txt

curl -skL "$1" |
  tee >(c/getURLs.js "$1" | grep -vxf ${OUT}/visited.txt >>${OUT}/urls.txt) |
  c/getText.js
