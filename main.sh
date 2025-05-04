#!/bin/bash

error() {
    echo "Error: $1" >/dev/stderr
    exit 1
}

correct() { [ "$(cut -d' ' -f 2 <"$BENCHMARK.hash" | grep -c 1)" -eq 0 ]; }

is_integer() { [[ $1 =~ ^[0-9]+$ && $1 -gt 0 ]]; }

main() {
    if [[ $# -lt 1 ]]; then
        error "Usage: $0 BENCHMARK_NAME [--time|--resources|--bare|args...]"
    fi

    BENCHMARK=$(basename "$1")
    shift
    export BENCHMARK

    export LC_ALL=C
    KOALA_SHELL=${KOALA_SHELL:-bash}
    export KOALA_SHELL
    measure_time=false
    measure_resources=false
    run_locally=false
    runs=1
    stats_prefix="${KOALA_SHELL%% *}_${BENCHMARK}_stats"
    time_values=()
    stats_files=()
    args=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
        --time | -t)
            measure_time=true
            shift
            ;;
        --resources)
            measure_resources=true
            shift
            ;;
        --bare)
            run_locally=true
            shift
            ;;
        --runs | -n)
            shift
            [[ $# -eq 0 ]] && error "Missing value for -n/--runs"
            is_integer "$1" || error "Value for -n/--runs must be a positive integer"
            runs="$1"
            shift
            ;;
        *)
            args+=("$1")
            shift
            ;;
        esac
    done

    REPO_TOP=$(git rev-parse --show-toplevel)

    cd "$(dirname "$0")/$BENCHMARK" || error "Could not cd into benchmark folder"

    for ((i = 1; i <= runs; i++)); do
        # Download dependencies
        if ((i == 1)); then
            ./install.sh "${args[@]}" ||
                error "Failed to download dependencies for $BENCHMARK"
        fi
        # Fetch inputs
        if ((i == 1)) || [[ "$BENCHMARK" == "riker" ]]; then
            ./fetch.sh "${args[@]}" ||
                error "Failed to fetch inputs for $BENCHMARK"
        fi

        if [[ "$measure_resources" == true && "$run_locally" == false ]]; then
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
            mv -f "$REPO_TOP/$BENCHMARK/benchmark_stats.txt" \
                "$REPO_TOP/$BENCHMARK/${stats_prefix}.txt"

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

            mkdir -p logs
            pidstat_log="logs/${BENCHMARK}.pidstat"
            io_log="logs/${BENCHMARK}.io"
            mem_log="logs/${BENCHMARK}.mem"
            cpu_log="logs/${BENCHMARK}.cpu"
            time_log="logs/${BENCHMARK}.time"

            # Start monitoring tools
            (while true; do
                cat /proc/diskstats >>"$io_log"
                sleep "$interval"
            done) &
            IO_PID=$!

            (pidstat -rud -h "$interval" >"$pidstat_log" 2>/dev/null) &
            PIDSTAT_PID=$!

            (while true; do
                cat /proc/stat >>"$cpu_log"
                cat /proc/meminfo >>"$mem_log"
                sleep "$interval"
            done) &
            BACKUP_PID=$!

            trap 'kill -TERM $IO_PID $PIDSTAT_PID $BACKUP_PID 2>/dev/null' EXIT

            time_fmt=$'USER=%U\nSYS=%S\nELAPSED=%e\nMAXRSS=%M\nIOIN=%I\nIOOUT=%O'
            /usr/bin/time -f "$time_fmt" -o "$time_log" \
                ./execute.sh "${args[@]}" \
                1>"${BENCHMARK}.out" \
                2>"${BENCHMARK}.err"
            CMD_STATUS=$?

            kill -TERM $IO_PID $PIDSTAT_PID $BACKUP_PID 2>/dev/null || true
            wait $IO_PID $PIDSTAT_PID $BACKUP_PID 2>/dev/null || true

            while IFS='=' read -r k v; do declare "$k=$v"; done <"$time_log"

            cpu_total=$(awk -v u="$USER" -v s="$SYS" 'BEGIN{printf "%.2f", u+s}')
            io_bytes=$(awk -v io_in="$IOIN" -v io_out="$IOOUT" 'BEGIN{print (io_in + io_out) * 1024}')
            mem_bytes=$((MAXRSS * 1024))

            get_size() {
                local p=$1
                if [[ -f $p ]]; then
                    stat -c%s -- "$p"
                elif [[ -d $p ]]; then
                    du -sb -- "$p" | awk '{print $1}'
                else
                    return 1
                fi
            }
            if [[ -n $BENCHMARK_INPUT_FILE ]]; then
                if INPUT_BYTES=$(get_size "$BENCHMARK_INPUT_FILE"); then
                    : # success
                else
                    echo "Warning: $BENCHMARK_INPUT_FILE not found; byte metrics = 0" >&2
                    INPUT_BYTES=0
                fi
            else
                INPUT_BYTES=0
            fi

            cpu_per_byte=$(awk -v c="$cpu_total" -v b="$INPUT_BYTES" \
                'BEGIN{printf "%.6f", (b?c/b:0)}')
            mem_per_byte=$(awk -v m="$mem_bytes" -v b="$INPUT_BYTES" \
                'BEGIN{printf "%.6f", (b?m/b:0)}')
            io_per_byte=$(awk -v i="$io_bytes" -v b="$INPUT_BYTES" \
                'BEGIN{printf "%.6f", (b?i/b:0)}')

            stats_file="${stats_prefix}.txt"
            cat >"$stats_file" <<EOF
Benchmark Statistics
==================================================

Benchmark: $BENCHMARK
--------------------------------------------------
Total CPU time:        ${cpu_total} sec
Total Wall time:       ${ELAPSED} sec
Total IO bytes:        ${io_bytes}
Max Memory Usage:      ${mem_bytes} bytes
Total Input bytes:     ${INPUT_BYTES}
CPU time per input byte: ${cpu_per_byte} sec/byte
Memory per input byte:  ${mem_per_byte} bytes/byte
IO per input byte:      ${io_per_byte} bytes/byte
Time in Shell: 0.00 sec (not measured with --bare flag)
Time in Commands: ${cpu_total} sec
==================================================
EOF
            echo "Saved local stats to $stats_file"

        elif $measure_time; then
            if ! command -v /usr/bin/time &>/dev/null; then
                echo "Installing /usr/bin/time and gawk..."
                sudo apt-get update && sudo apt-get install -y time gawk
            fi

            echo "Timing benchmark: $BENCHMARK  (run #$i)"

            time_val_file="${BENCHMARK}_time_run${i}.txt"
            rm -f "$time_val_file"

            /usr/bin/time -f "%e" -o "$time_val_file" \
                ./execute.sh "${args[@]}" \
                1>"${BENCHMARK}.out" \
                2>"${BENCHMARK}.err"
            CMD_STATUS=$?

            if [[ -s "$time_val_file" ]]; then
                runtime=$(<"$time_val_file")
            else
                echo "Warning: could not capture runtime for run #$i" >&2
                runtime=0
            fi

            time_values+=("$runtime")
            [[ $CMD_STATUS -ne 0 ]] && error "Failed to run $BENCHMARK"

        else
            ./execute.sh "${args[@]}" >"$BENCHMARK.out" 2>"$BENCHMARK.err" || error "Failed to run $BENCHMARK"
        fi

        # Verify output
        ./validate.sh "${args[@]}" >"$BENCHMARK.hash" || error "Failed to verify output for $BENCHMARK"

        # Cleanup outputs
        if (( i == runs )) || [[ "$BENCHMARK" == "riker" ]]; then
            ./clean.sh "${args[@]}"
        fi

        if correct; then
            echo "$BENCHMARK [pass]"
        else
            error "$BENCHMARK [fail]"
        fi

        if [[ $measure_resources == true ]]; then
            src_stats="$REPO_TOP/$BENCHMARK/${stats_prefix}.txt"

            if [[ -f $src_stats ]]; then
                dst_stats="$REPO_TOP/$BENCHMARK/${stats_prefix}_run${i}.txt"
                cp -f -- "$src_stats" "$dst_stats"
                stats_files+=("$dst_stats") # remember it for later aggregation
            else
                echo "Warning: $src_stats not found for run #$i" >&2
            fi
        fi
    done

    if [[ $measure_resources == true && "$run_locally" == false && ${#stats_files[@]} -gt 1 ]]; then
        agg_script="$REPO_TOP/infrastructure/aggregate_stats.py"
        if [[ -f $agg_script ]]; then
            python3 "$agg_script" "${stats_files[@]}" \
                >"$REPO_TOP/$BENCHMARK/${stats_prefix}_aggregated.txt" ||
                echo "Aggregation failed" >&2
            echo "Wrote aggregated stats to ${stats_prefix}_aggregated.txt"
        else
            echo "Aggregation script $agg_script missing" >&2
        fi
    fi

    if $measure_time && ((${#time_values[@]} > 1)); then
        times_file="${BENCHMARK}_times_aggregated.txt"

        {
            echo "Aggregated Wall-Clock Runtimes"
            echo "========================================"
            printf "Runs analysed: %s\n\n" "${#time_values[@]}"

            printf "%s\n" "${time_values[@]}" |
                awk '
                { sum += $1; arr[NR] = $1 }
                END {
                    asort(arr);                       # gawk ≥ 4
                    mean = sum / NR
                    printf "Mean  : %.3f sec\n", mean
                    printf "Min   : %.3f sec\n", arr[1]
                    printf "Max   : %.3f sec\n", arr[NR]
                }
            '

            echo
            echo "Per-run raw values:"
            paste <(seq 1 ${#time_values[@]}) <(printf "%s\n" "${time_values[@]}") |
                awk '{printf "  run %-3s : %.3f sec\n", $1, $2}'
        } >"$times_file"

        echo "Wrote aggregated runtimes to $times_file"
    fi

    cd - || exit 1

}

main "$@"
