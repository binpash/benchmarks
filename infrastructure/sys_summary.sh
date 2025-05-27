#!/bin/bash

output_file="benchmark_results.csv"

default_benchmarks=(
  "aurpkg"
  "bio"  
  "covid-mts"
  "dpt"
  "file-enc"
  "git-workflow"
  "llm"
  "log-analysis"
  "makeself"
  "max-temp"
  "media-conv"
  "nlp"
  "oneliners"
  "prog-inf"
  "riker"
  "sklearn"
  "teraseq"
  "tuft-weather"
  "unix50"
  "vps-audit"
  "web-index"
)

benchmarks=()
args=()

for arg in "$@"; do
    if [[ "$arg" == --* ]]; then
        args+=("$arg")
    else
        benchmarks+=("$arg")
    fi
done

if [ ${#benchmarks[@]} -eq 0 ]; then
    benchmarks=("${default_benchmarks[@]}")
fi

error() {
    echo "Error: $1" >&2
    exit 1
}

echo "Benchmark,Sys Calls,File Descriptors" > "$output_file"

for benchmark in "${benchmarks[@]}"; do
    echo "Running benchmark: $benchmark"

    mkdir -p "/tmp/${benchmark}"
    strace_output="/tmp/${benchmark}_strace.txt"
    lsof_output="/tmp/${benchmark}_lsof.txt"

    if ! cd "./$benchmark"; then
        echo "$benchmark,FAIL: directory not found,0" >> "../$output_file"
        continue
    fi

    setsid ./execute.sh "${args[@]}" &
    pid=$!

    pgid=$(ps -o pgid= -p "$pid" | tr -d ' ')
    sleep 1

    if [[ -n "$pgid" ]]; then
        lsof -g "$pgid" > "$lsof_output" || echo "Warning: lsof failed"
    else
        echo "Warning: Failed to get PGID for $pid"
        > "$lsof_output"
    fi

    wait "$pid"

    strace -c -f -o "$strace_output" ./execute.sh "${args[@]}" || {
        echo "$benchmark,FAIL: strace failed,0" >> "../$output_file"
        cd - > /dev/null || exit
        continue
    }

    total_syscalls=$(awk '/^100.00/ {print $4}' "$strace_output")
    total_syscalls=${total_syscalls:-0}

    fd_count=$(wc -l < "$lsof_output")
    fd_count=${fd_count:-0}

    echo "$benchmark,$total_syscalls,$fd_count" >> "../$output_file"

    rm -f "$strace_output" "$lsof_output"
    cd - > /dev/null || exit
done

echo "Benchmark results saved to $output_file."
