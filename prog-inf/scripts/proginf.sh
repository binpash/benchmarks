#!/bin/bash
TOP=$(realpath "$(dirname "$0")/..")
MIR_BIN=$TOP/inputs/mir-sa/.bin/mir-sa
OUT=$TOP/outputs
IN=$TOP/inputs
INDEX=${INDEX:-"$TOP/inputs/index.txt"}

mkdir -p "${OUT}/"
pkg_count=0
while read -r package
do
    pkg_count=$((pkg_count + 1));
    cd "$IN/node_modules/$package" || exit 1
    ${MIR_BIN} -p  > "${OUT}/$pkg_count.log"
done < "$INDEX"
