#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-pcap"
results_dir="${eval_dir}/results"
hashes_dir="${eval_dir}/hashes"

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd

if [[ "$@" == *"--generate"* ]]; then
    md5sum pcap_bench$suffix/* > $hashes_dir/pcap_bench$suffix.md5sum
    md5sum count_packets$suffix/* > $hashes_dir/count_packets$suffix.md5sum
fi

okay=0
if ! md5sum --check --quiet $hashes_dir/pcap_bench$suffix.md5sum; then
    okay=1
    echo "pcap_bench $suffix failed verification"
fi
if ! md5sum --check --quiet $hashes_dir/count_packets$suffix.md5sum; then
    okay=1
    echo "count_packets $suffix failed verification"
fi

exit $okay
