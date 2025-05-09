#!/bin/bash

export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

output_dir="$PASH_SPEC_TOP/report/output/max_temp/small" # Adjust the path as necessary
download_dir="$PASH_SPEC_TOP/report/resources/max_temp"
result_dir="$PASH_SPEC_TOP/results/max_temp/small"
echo $output_dir
echo $download_dir
echo $result_dir

rm -f $download_dir/*.txt

./run --target sh-only
rm $download_dir/*.txt
mkdir -p $result_dir/base/
mv $output_dir/* $result_dir/base/

for window in 40 30 20 10 0
do
    ./run --target hs-only --window $window
    rm $download_dir/*.txt
    mkdir -p $result_dir/$window
    mv $output_dir/* $result_dir/$window
done

