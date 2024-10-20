REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-pcap"
results_dir="${eval_dir}/results"
input_dir="${eval_dir}/input"
scripts_dir="${eval_dir}/scripts"

mkdir -p $results_dir

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

export INPUT=${input_dir}/201011271400.dump
export INPUT2=${input_dir}/2018-07-20-17-31-20-192.168.100.108.pcap

script="${scripts_dir}/pcap_bench.sh"
out_dir=$results_dir/pcap_bench$suffix
mkdir -p $out_dir
$script > $out_dir/out

script="${scripts_dir}/count_packets.sh"
out_dir=$results_dir/count_packets$suffix
mkdir -p $out_dir
$script > $out_dir/out
