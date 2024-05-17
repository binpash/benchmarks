REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/analysis-pcap"
results_dir="${eval_dir}/results"
inputs_dir="${eval_dir}/input"

shell="/bin/bash"
mkdir -p $results_dir

# script2="${eval_dir}/split_pcap.sh"
# script3="${eval_dir}/count_packets.sh"

export IN1=${inputs_dir}/wav
export IN2=${inputs_dir}/rtf
export IN3=${inputs_dir}/jpg
export IN4=${inputs_dir}/linux
export IN5=${inputs_dir}/
export IN7=${inputs_dir}/

export OUT1=${results_dir}/out1
export OUT3=${results_dir}/out3
export OUT6=${results_dir}/out6

# export OUTPUT=

script="${eval_dir}/pcap_bench.sh"
echo "Executing $(basename "$script")"
"$shell" "$script" > "$results_dir/$(basename "$script").out"

script="${eval_dir}/count_packets.sh"
echo "Executing $(basename "$script")"
"$shell" "$script" > "$results_dir/$(basename "$script").out"
