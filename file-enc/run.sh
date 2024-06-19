# Run benchmarks
# ~ Execute & hash outputs
# ~ Add metrics of performance
eval_dir="./scripts"
in_dir="./inputs/pcap_data"
out_dir="./outputs"
shell="/bin/bash"

if [[ "$1" == "--small" ]]; then
    echo "Using small input"
    in_dir="./inputs/pcap_data"
else
    echo "Using default input"
    in_dir="./inputs/pcap_data_small"
fi


scripts=(
    "compress_files"
    "encrypt_files"
)

file-enc() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/file-enc.res"
    > $mode_res_file
}

cd "$(realpath $(dirname "$0"))"
mkdir -p $out_dir

script="${eval_dir}/compress_files.sh"
echo "Executing $(basename "$script")"
scriptname=$(basename "script")
time_file="$out_dir/$scriptname.time"
time "$shell" "$script" "$in_dir" "$out_dir/compress" > "$out_dir/$(basename "$script").out" 2> $time_file

# script="${eval_dir}/dependencies.sh"
# echo "Executing $(basename "$script")"
# "$shell" "$script" > "$out_dir/$(basename "$script").out"