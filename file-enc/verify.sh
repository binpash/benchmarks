# The output files .enc and .zip will vary across separate runs
# even for the same input file. Therefore, correctness check is omitted.

#md5sum compress_files.full/* > /benchmarks/file-enc/hashes/compress_files.full.md5sum
#md5sum encrypt_files.full/* > /benchmarks/file-enc/hashes/encrypt_files.full.md5sum
#md5sum compress_files.full/* > /benchmarks/file-enc/hashes/compress_files.small.md5sum
#md5sum encrypt_files.full/* > /benchmarks/file-enc/hashes/encrypt_files.small.md5sum

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
results_dir="${eval_dir}/results"
hashes_dir="${eval_dir}/hashes"

suffix=".full"
if [[ $1 == "--small" ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd
okay=0
if ! md5sum --check --quiet $hashes_dir/encrypt_files$suffix.md5sum; then
    okay=1
    echo "encrypt_files $suffix failed verification"
fi
if ! md5sum --check --quiet $hashes_dir/compress_files$suffix.md5sum; then
    okay=1
    echo "compress_files $suffix failed verification"
fi

exit $okay
