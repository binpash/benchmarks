REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-pcap"
results_dir="${eval_dir}/results"
inputs_dir="${eval_dir}/input"

shell="/bin/bash"
mkdir -p $results_dir

# script2="${eval_dir}/split_pcap.sh"
# script3="${eval_dir}/count_packets.sh"

export INPUT=${inputs_dir}/201011271400.dump
export INPUT2=${inputs_dir}/2018-07-20-17-31-20-192.168.100.108.pcap

# export OUTPUT=

script="${eval_dir}/pcap_bench.sh"
echo "Executing $(basename "$script")"
"$shell" "$script" > "$results_dir/$(basename "$script").out"

script="${eval_dir}/count_packets.sh"
echo "Executing $(basename "$script")"
"$shell" "$script" > "$results_dir/$(basename "$script").out"
