#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)

tests="5TERA 5TERA-short 5TERA3 TERA3 Akron5Seq mouse_SIRV dRNASeq RNASeq RiboSeq"

for test in ${tests}
do
	(testdir="${REPO_ROOT}/teraseq/${test}"; echo "downloading for ${test}"; cd ${testdir}; bash download.sh) &
done
