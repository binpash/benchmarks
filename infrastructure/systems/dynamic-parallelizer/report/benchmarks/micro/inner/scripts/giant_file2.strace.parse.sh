#! /bin/sh
generate_unique_file() {
    local dir="$OUTPUT_DIR"
    local prefix="strace_log"
    local counter_file="$dir/${prefix}"
    if [ ! -f "$counter_file" ]; then
        echo 0 > "$counter_file"
    fi
    local counter
    counter=$(<"$counter_file")
    counter=$((counter + 1))
    echo "$counter" > "$counter_file"
    local filename="$dir/${prefix}_${counter}"
    echo "$filename"
} 
hs_base=$(git rev-parse --show-toplevel)
PARSE="python3 $hs_base/parallel-orch/trace_v2.py"


OUTPUT=${OUTPUT:-.}
hs_base=$(git rev-parse --show-toplevel)
SCRIPTS="${hs_base}/report/benchmarks/micro/scripts"
touch "$OUTPUT"/giant
logfile=$(generate_unique_file)
strace -y -f  --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i python3 "$SCRIPTS"/giant_file.py "$OUTPUT"/giant 10000
$PARSE $logfile > $(generate_unique_file)
