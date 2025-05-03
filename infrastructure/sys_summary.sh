#!/bin/bash

output_file="benchmark_results.csv"

benchmarks=(
  "aurpkg"
  "bio"
  "covid-mts"
  "file-enc"
  "log-analysis"
  "makeself"
  "max-temp"
  "media-conv"
  "nlp"
  "oneliners"
  "riker"
  "sklearn"
  "unix50"
  "vps-audit"
  "web-index"
)

# Error handler
error() {
    echo "Error: $1" >&2
    exit 1
}

echo "Benchmark,Sys Calls,File Descriptors" > "$output_file"

for benchmark in "${benchmarks[@]}"; do
    echo "Running benchmark: $benchmark"

    strace_output="/tmp/${benchmark}_strace.txt"
    lsof_output="/tmp/${benchmark}_lsof.txt"

    if ! cd "./$benchmark"; then
        echo "$benchmark [fail]: Directory not found" >> "$output_file"
        continue
    fi

    strace -c -o "$strace_output" ./execute.sh --small || {
        echo "$benchmark [fail]: strace failed" >> "$output_file"
        cd - > /dev/null
        continue
    }

    ./execute.sh --small &
    pid=$!

    sleep 1

    if [[ -d /proc/$pid ]]; then
        lsof -p "$pid" > "$lsof_output"
    else
        echo "Warning: Process $pid ended before lsof could capture file descriptors."
        > "$lsof_output"
    fi

    wait "$pid"

    total_syscalls=$(awk '/^100.00/ {print $4}' "$strace_output")
    total_syscalls=${total_syscalls:-0}

    fd_count=$(wc -l < "$lsof_output")
    fd_count=${fd_count:-0}

    echo "$benchmark,$total_syscalls,$fd_count" >> "../$output_file"

    rm -f "$strace_output" "$lsof_output"
    cd - > /dev/null
done

echo "Benchmark results saved to $output_file."
