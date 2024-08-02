#!/bin/bash
# encrypt all files in a directory 
mkdir -p $2

pure_func() {
    openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -k 'key'
}
export -f pure_func

# item=$1
# output_name=$(basename $item).enc
# hdfs dfs -cat -ignoreCrc $item | pure_func > $2/$output_name

for item in $(ls ${1});
do
    output_name=$(basename $item).enc
    cat $item | pure_func > $2/$output_name
done

echo 'done';
