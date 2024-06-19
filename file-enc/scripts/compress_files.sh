#!/bin/bash
# compress all files in a directory
IN=$1
OUT=$2
LOGS=${OUT}/logs
mkdir -p $LOGS
run_tests() {
    name=$(basename $1).zip
    zip -r ${OUT}/$name $1
}

export -f run_tests

pkg_count=0
for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > "${LOGS}"/"$pkg_count.log"
done

echo 'done';