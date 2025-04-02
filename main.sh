#!/bin/bash

error() {
    echo "Error: $1" >/dev/stderr
    exit 1
}

correct() { [ "$(cut -d' ' -f 2 < "$BENCHMARK.hash" | grep -c 1)" -eq 0 ]; }

main() {
    if [[ $# -lt 1 ]]; then
        error "Usage: $0 BENCHMARK_NAME [--time|--resources|--bare|args...]"
    fi

    BENCHMARK=$(basename "$1")
    shift
    export BENCHMARK

    measure_time=false
    measure_resources=false
    run_locally=false

    args=()
    for arg in "$@"; do
        case "$arg" in
        --time | -t)
            measure_time=true
            ;;
        --resources)
            measure_resources=true
            ;;
        --bare)
            run_locally=true
            ;;
        *)
            args+=("$arg")
            ;;
        esac
    done

    REPO_TOP="/benchmarks"

    cd "$(dirname "$0")/$BENCHMARK" || error "Could not cd into benchmark folder"

    # Download dependencies
    ./deps.sh "${args[@]}" || error "Failed to download dependencies for $BENCHMARK"

    # Fetch inputs
    ./input.sh "${args[@]}" || error "Failed to fetch inputs for $BENCHMARK"

    if [[ "$measure_resources" == true && "$run_locally" == false ]]; then
        #! for some reason some benchmark results are being tainted when they are ran through the dynamic analysis
        #! resourse stats look good though
        #!TODO check why this is happening

        echo "[*] Running dynamic resource analysis for $BENCHMARK"
        sudo apt-get install -y autoconf automake libtool build-essential cloc
        pip install --break-system-packages -r "$REPO_TOP/infrastructure/requirements.txt"

        mkdir -p "$REPO_TOP/infrastructure/target/process-logs"
        mkdir -p "$REPO_TOP/infrastructure/target/backup-process-logs"
        find "$REPO_TOP/infrastructure/target/process-logs" -type f \
            -exec mv {} "$REPO_TOP/infrastructure/target/backup-process-logs/" \; || true
        rm -f "$REPO_TOP"/infrastructure/target/process-logs/*
        rm -f "$REPO_TOP"/infrastructure/target/dynamic_analysis.jsonl

        cd "$REPO_TOP" || exit 1
        python3 "$REPO_TOP/infrastructure/run_dynamic.py" "$BENCHMARK" "${args[@]}" || error "Failed to run $BENCHMARK"

        cd "$REPO_TOP/infrastructure" || exit 1
        make target/dynamic_analysis.jsonl
        python3 viz/dynamic.py "$REPO_TOP/$BENCHMARK" --text
        cat "$REPO_TOP/$BENCHMARK/benchmark_stats.txt"

        cd "$REPO_TOP/$BENCHMARK" || exit 1

    elif [[ "$measure_resources" == true && "$run_locally" == true ]]; then
        echo "Running local resource monitoring for $BENCHMARK"
        interval=0.1

        if ! command -v pidstat &>/dev/null; then
            echo "Installing pidstat..."
            sudo apt-get update && sudo apt-get install -y sysstat
        fi

        if ! command -v /usr/bin/time &>/dev/null; then
            echo "Installing /usr/bin/time..."
            sudo apt-get update && sudo apt-get install -y time
        fi

        pidstat_log="$BENCHMARK.pidstat"
        io_log="$BENCHMARK.io"
        mem_log="$BENCHMARK.mem"
        cpu_log="$BENCHMARK.cpu"
        mkdir -p logs

        # Start monitoring tools
        (while true; do
            cat /proc/diskstats >>"logs/$io_log"
            sleep "$interval"
        done) &
        IO_PID=$!

        (pidstat -rud -h "$interval" >"logs/$pidstat_log" 2>/dev/null) &
        PIDSTAT_PID=$!

        (while true; do
            cat /proc/stat >>"logs/$cpu_log"
            cat /proc/meminfo >>"logs/$mem_log"
            sleep "$interval"
        done) &
        BACKUP_PID=$!

        trap 'kill -TERM $IO_PID $PIDSTAT_PID $BACKUP_PID 2>/dev/null' EXIT

        # Run benchmark and measure time in background
        /usr/bin/time -v ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" &
        CMD_PID=$!
        wait "$CMD_PID"


        # Stop monitoring tools
        kill -TERM $IO_PID $PIDSTAT_PID $BACKUP_PID 2>/dev/null || true
        wait $IO_PID $PIDSTAT_PID $BACKUP_PID 2>/dev/null || true

    elif $measure_time; then
        echo "Timing benchmark: $BENCHMARK"
        /usr/bin/time -f "Runtime: %E (CPU: %P)" ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
    else
        ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
    fi

    # Verify output
    ./verify.sh "${args[@]}" >"$BENCHMARK.hash" || error "Failed to verify output for $BENCHMARK"

    # Cleanup outputs
    ./cleanup.sh "${args[@]}"

    if correct; then
        echo "$BENCHMARK [pass]"
    else
        error "$BENCHMARK [fail]"
    fi

    cd - || exit 1
}

main "$@"
