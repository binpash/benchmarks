#!/bin/bash

error() {
    echo "Error: $1" >/dev/stderr
    exit 1
}

correct() { [ "$(cut -d' ' -f 2 <"$BENCHMARK.hash" | grep -vc 0)" -eq 0 ]; }

in_container() {
  [ -f /proc/1/cgroup ] && grep -qaE 'docker|kubepods|containerd' /proc/1/cgroup && return 0
  [ -f /.dockerenv ] && return 0
  return 1
}

is_integer() { [[ $1 =~ ^[0-9]+$ && $1 -gt 0 ]]; }

usage() {
    echo "Usage: $0 BENCHMARK_NAME [--time|--resources|--bare|args...]"
    echo "  --min            Run the benchmark with minimal inputs (default)"
    echo "  --small          Run the benchmark with small inputs"
    echo "  --full          Run the benchmark with full inputs"
    echo "  --time, -t       Measure wall-clock time"
    echo "  --resources      Measure resource usage"
    echo "  --bare           Run locally without Docker"
    echo "  --runs, -n N     Number of runs (default: 1)"
    echo "  --clean, -c      Run the full cleanup script (both inputs and outputs)"
    echo "  --keep, -k       Keep outputs"
    echo "  --prune          Run the benchmark on a fresh container (will need to re-download everything on each run)"
    echo "  --help, -h       Show this help message"
}

main() {
    if [[ $# -lt 1 ]]; then
        usage
        exit 1
    fi


    measure_time=false
    measure_resources=false
    run_locally=false
    run_cleanup=false
    keep_outputs=false
    prune=false
    runs=1
    size="min"

    args=()
    main_args=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
        --help | -h) 
            usage
            exit 0
            ;;
        --time | -t)
            measure_time=true
            main_args+=("$1")
            shift
            ;;
        --resources)
            measure_resources=true
            main_args+=("$1")
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
        --clean | -c)
            run_cleanup=true
            main_args+=("$1")
            shift
            ;;
        --keep | -k)
            keep_outputs=true
            main_args+=("$1")
            shift
            ;;
        --prune)
            prune=true
            main_args+=("$1")
            shift
            ;;
        --min)
            size="min"
            main_args+=("$1")
            shift
            ;;
        --small)
            size="small"
            main_args+=("$1")
            shift
            ;;
        --full)
            size="full"
            main_args+=("$1")
            shift
            ;;
        *)
            if [[ "$1" != -* ]]; then
                BENCHMARK="$(basename "$1")"
            else
                args+=("$1")
            fi
            shift
            ;;
        esac
    done

    if [[ "$size" == "min" ]]; then
        args+=("--min")
    elif [[ "$size" == "small" ]]; then
        args+=("--small")
    elif [[ "$size" == "full" ]]; then
        args+=("")
    fi

    [ -z "$BENCHMARK" ] && usage && exit 1
    [ ! -d "$BENCHMARK" ] && error "Benchmark folder $BENCHMARK does not exist"
    export BENCHMARK

    export LC_ALL=C
    export TZ=UTC
    KOALA_SHELL=${KOALA_SHELL:-bash}
    KOALA_CONTAINER_CMD=${KOALA_CONTAINER_CMD:-docker}
    export KOALA_SHELL
    export KOALA_INFO="time:mem:io:cpu"

    if  in_container; then
        run_locally=true
    fi

    shell_word=${KOALA_SHELL%% *}
    shell_word=${shell_word##*/}
    shell_safe=${shell_word//[^A-Za-z0-9_.-]/_}
    echo "Using shell: $KOALA_SHELL"
    stats_prefix="${BENCHMARK}_${shell_safe}_stats"
    time_values=()
    stats_files=()

    REPO_TOP=$(git rev-parse --show-toplevel)
    [[ -z "$REPO_TOP" ]] && error "Failed to determine repository top"
    
    if ! $run_locally; then
        DOCKER_IMAGE=${KOALA_DOCKER_IMAGE:-ghcr.io/binpash/benchmarks:latest}
        echo "Launching KOALA in $KOALA_CONTAINER_CMD container ($DOCKER_IMAGE)"
        $KOALA_CONTAINER_CMD pull "$DOCKER_IMAGE"

        if $prune; then
            echo "Running with prune mode: starting clean container"
            $KOALA_CONTAINER_CMD run --rm \
                "$DOCKER_IMAGE" \
                -w "/benchmarks" \
                ./main.sh "$BENCHMARK" ${args[*]} ${main_args[*]} --bare
        else
            echo "Mounting $REPO_TOP to /benchmarks in the container"
            $KOALA_CONTAINER_CMD run --rm \
                -v "$REPO_TOP":/benchmarks \
                -w "/benchmarks" \
                -e KOALA_SHELL="$KOALA_SHELL" \
                "$DOCKER_IMAGE" \
                ./main.sh "$BENCHMARK" ${args[*]} ${main_args[*]} --bare
        fi
        exit $?
    fi

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

        # Delete outputs before each run
        if [ "$BENCHMARK" != "riker" ]; then
            ./clean.sh "${args[@]}"
        fi

        echo "Executing $BENCHMARK $(date) ($i/$runs)"
        if [[ "$measure_resources" == true ]]; then
            echo "[*] Running dynamic resource analysis for $BENCHMARK"
            # check if deps are installed
            if ! command -v cloc &>/dev/null || ! command -v python3 &>/dev/null; then
                echo "Please run setup.sh first to install dependencies."
                exit 1
            fi

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
            python3 viz/dynamic.py "$REPO_TOP/$BENCHMARK" # --text
            if [[ -f "$REPO_TOP/$BENCHMARK/benchmark_stats.txt" ]]; then
                cat "$REPO_TOP/$BENCHMARK/benchmark_stats.txt"
                mv -f "$REPO_TOP/$BENCHMARK/benchmark_stats.txt" \
                "$REPO_TOP/$BENCHMARK/${stats_prefix}.txt" || error "Failed to move benchmark stats"
            else
                error "Failed to generate benchmark stats"
            fi

            find "$REPO_TOP/infrastructure/target/backup-process-logs" -type f \
                -exec mv {} "$REPO_TOP/infrastructure/target/process-logs/" \; || true
            cd "$REPO_TOP/$BENCHMARK" || exit 1

        elif $measure_time; then
            if ! command -v /usr/bin/time &>/dev/null || ! command -v gawk &>/dev/null; then
                echo "Please run setup.sh first to install dependencies."
                exit 1
            fi

            echo "Timing benchmark: $BENCHMARK  (run #$i)"

            time_val_file="${BENCHMARK}_${shell_safe}_time_run${i}.txt"
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
            ./execute.sh "${args[@]}" 2>"$BENCHMARK.err" | tee "$BENCHMARK.out" || error "Failed to run $BENCHMARK"
        fi

        # Verify output
        ./validate.sh "${args[@]}" >"$BENCHMARK.hash" || error "Failed to verify output for $BENCHMARK"

        # Cleanup outputs
        if [ "$keep_outputs" = false ] && [ "$i" -eq "$runs" ]; then
            if [ "$run_cleanup" = true ]; then
                ./clean.sh -f "${args[@]}"
            else
                ./clean.sh "${args[@]}"
            fi
        fi

        if correct; then
            echo "$BENCHMARK [pass]"
        else
            echo "$BENCHMARK [fail]"
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

cd "$(dirname "$0")" || error "Could not cd into script folder"

main "$@"
