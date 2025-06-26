#!/bin/bash

SAMPLING_INT=1
CSV_HEADER="Benchmark,Sys Calls,FD_Snapshot,Unique_FD,Peak_FD"

output_file="benchmark_results.csv"

default_benchmarks=(
  "analytics"
  "bio"  
  "ci-cd"
  "covid"
  "file-mod"
  "inference"
  "ml"
  "nlp"
  "oneliners"
  "pkg"
  "repl"
  "unixfun"
  "weather"
  "web-search"
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

[[ ${#benchmarks[@]} -eq 0 ]] && benchmarks=("${default_benchmarks[@]}")

echo "$CSV_HEADER" > "$output_file"

for benchmark in "${benchmarks[@]}"; do
    echo "Running benchmark: $benchmark"

    mkdir -p "/tmp/${benchmark}"
    strace_out="/tmp/${benchmark}_strace.txt"
    snap_out="/tmp/${benchmark}_lsof_snapshot.txt"
    stream_out="/tmp/${benchmark}_lsof_stream.txt"

    if ! cd "./$benchmark"; then
        echo "$benchmark,FAIL: directory not found" >> "../$output_file"
        continue
    fi

    ./install.sh || { echo "$benchmark,FAIL: install.sh" >> "../$output_file"; cd - >/dev/null; continue; }
    ./fetch.sh "${args[@]}" || { echo "$benchmark,FAIL: fetch.sh" >> "../$output_file"; cd - >/dev/null; continue; }

    setsid ./execute.sh "${args[@]}" &
    pid=$!
    sleep 1
    pgid=$(ps -o pgid= -p "$pid" | tr -d ' ')

    [[ -z $pgid ]] && { echo "$benchmark,FAIL: no PGID" >> "../$output_file"; cd - >/dev/null; continue; }

    lsof -n -P -w -g "$pgid" -r${SAMPLING_INT} > "$stream_out" &
    sampler_pid=$!

    lsof -n -P -w -g "$pgid" > "$snap_out" || echo "Warning: lsof snapshot failed"

    wait "$pid"
    kill "$sampler_pid" 2>/dev/null

    strace -c -f -o "$strace_out" ./execute.sh "${args[@]}" || {
        echo "$benchmark,FAIL: strace" >> "../$output_file"; cd - >/dev/null; continue; }

    syscalls=$(awk '/^100\.00/ {print $4}' "$strace_out")

    FD_Snapshot=$(awk 'NR>1' "$snap_out" | wc -l)
    Unique_FD=$(awk 'NR>1 {print $4":"$10}' "$snap_out" | sort -u | wc -l)

    if [[ -s "$stream_out" ]]; then
        Peak_FD=$(awk '
            NF && $1!="COMMAND" {seen[$4":"$10]++}
            /^====/             {print length(seen); delete seen}
        ' "$stream_out" | awk 'max<$1{max=$1} END{print max}')
    else
        Peak_FD=0
    fi

    echo "${benchmark%,/},${syscalls:-0},${FD_Snapshot:-0},${Unique_FD:-0},${Peak_FD:-0}" \
     >> "../$output_file"

    rm -f "$strace_out" "$snap_out" "$stream_out"
    cd - >/dev/null
done

echo "Benchmark results saved to $output_file."
