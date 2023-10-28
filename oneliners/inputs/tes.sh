#!/bin/bash    
touch "1M.txt"
while [ $(wc -c <"./1M.txt") -lt 1000000 ]
do
    cat ./temp/original.txt >> ./1M.txt 
done  