#!/bin/bash

while read line
do
pop="${s_line%%-*}"
temp="${s_line#*-}"
sample="${temp%%-*}"
ftp_url=$(echo $line | awk -F'-' '{print $NF}')
size_bytes=$(curl -sI $ftp_url | grep -i Content-Length | awk '{print $2}' | tr -d '\r')
size_gb=$(echo "scale=2; $size_bytes / 1000^3" | bc)
echo -n "$size_gb GB "
echo $line
done < input_all.txt
