#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/covid-mts"
results_dir="${eval_dir}/results"

md5hash=$(md5 -q $results_dir/{1,2,3,4,5}.out)
if [ "$md5hash" == "5e8e7590464ef2983333b421b03e427a
d99ee01e6a9caf54bf120ffc68e92ba1
7f7dd7e656255e21982660b2f765a5b0
ce97fb55b8b9703dbd7a968dcb26718f
d41d8cd98f00b204e9800998ecf8427e" ];
then
    echo "Valid"
else
    echo "Invalid"
    exit 1
fi
