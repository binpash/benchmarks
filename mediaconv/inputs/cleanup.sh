#!/bin/bash

cd $(dirname $0)

IN=${IN:./}
OUT=${OUT:./}

rm -rf ${IN}/jpg
rm -rf ${IN}/log_data
rm -rf ${IN}/wav
rm -rf ${IN}/nginx-logs
rm -rf ${IN}/node_modules
rm -rf ${IN}/pcap_data
rm -rf ${IN}/pcaps
rm -rf ${IN}/packages
rm -rf ${IN}/mir-sa
rm -rf ${IN}/deps
rm -rf ${IN}/bio
rm -rf ${IN}/output
rm -rf ${OUT}
exit
