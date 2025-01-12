#!/bin/bash

error() {
    echo "Error: $1" >&2
    exit 1
}

correct() { 
    [ "$(cat $BENCHMARK.hash | cut -d' ' -f 2 | grep -c 1)" -eq 0 ]
}

profile_run() {
    local benchmark_name="$1"
    local output_prefix="$2"

    echo "Profiling ./run.sh for $benchmark_name..."

    strace -c -o "${output_prefix}_strace.txt" ./run.sh --small || error "Failed to run $benchmark_name"

    ./run.sh --small &
    local pid=$!

    sleep 1

    if [[ -d /proc/$pid ]]; then
        lsof -p "$pid" > "${output_prefix}_lsof.txt"
    else
        echo "Warning: Process $pid ended before lsof could capture file descriptors."
        > "${output_prefix}_lsof.txt"
    fi

    wait "$pid"

    local total_syscalls
    total_syscalls=$(awk '/^100.00/ {print $4}' "${output_prefix}_strace.txt")
    echo "Total system calls for $benchmark_name: $total_syscalls"

    if [[ -s "${output_prefix}_lsof.txt" ]]; then
        local fd_count
        fd_count=$(wc -l < "${output_prefix}_lsof.txt")
        echo "Total open file descriptors for $benchmark_name: $fd_count"
    else
        echo "No file descriptors captured for $benchmark_name."
    fi

    echo "Finished profiling ./run.sh for $benchmark_name."
}


main() {
    export BENCHMARK="$1"
    shift

    cd "$(dirname "$0")/$BENCHMARK" || exit 1

    ./deps.sh "$@" || error "Failed to download dependencies for $BENCHMARK"
    ./input.sh "$@" || error "Failed to fetch inputs for $BENCHMARK"

    OUTPUT_DIR="./profiling_results"
    mkdir -p "$OUTPUT_DIR"
    local output_prefix="${OUTPUT_DIR}/${BENCHMARK}"

    profile_run "$BENCHMARK" "$output_prefix"

    ./verify.sh "$@" > "$BENCHMARK.hash" || error "Failed to verify output for $BENCHMARK"

    ./cleanup.sh "$@"

    if correct; then
        echo "$BENCHMARK [pass]"
    else
        error "$BENCHMARK [fail]"
    fi

    cd - > /dev/null || exit 1
}

main "$@"