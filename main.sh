#!/bin/bash

error()
{
    echo "Error: $1" > /dev/stderr
    exit 1
}

correct() { [ "$(cat $BENCHMARK.hash | cut -d' ' -f 2 | grep -c 1)" -eq 0 ]; }

main()
{
    export BENCHMARK="$1"
    shift

    measure_time=false
    measure_resources=false

    args=()
    for arg in "$@"; do
        case "$arg" in
            --time)
                measure_time=true
                ;;
            -t)
                measure_time=true
                ;;
            --resources)
                measure_resources=true
                ;;
            *)
                args+=("$arg")
                ;;
        esac
    done

    cd "$(dirname "$0")/$BENCHMARK" || exit 1

    # Download dependencies
    ./deps.sh "${args[@]}" || error "Failed to download dependencies for $BENCHMARK"

    # Fetch inputs
    ./input.sh "${args[@]}" || error "Failed to fetch inputs for $BENCHMARK"

    # Run benchmark with time/resource measurement
    if $measure_time && $measure_resources; then
        /usr/bin/time -f "Runtime: %E (CPU: %P)" ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" &
        time_pid=$!
        benchmark_pid=$!
        pidstat -d -r -u -p "$benchmark_pid" > "$BENCHMARK.resources"
        wait $time_pid
    elif $measure_time; then
        /usr/bin/time -f "Runtime: %E (CPU: %P)" ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
    elif $measure_resources; then
        ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
        benchmark_pid=$!
        pidstat -d -r -u -p "$benchmark_pid" > "$BENCHMARK.resources"

        wait $benchmark_pid
    else
        ./run.sh "${args[@]}" > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
    fi

    # Verify output
    ./verify.sh "${args[@]}" > "$BENCHMARK.hash" || error "Failed to verify output for $BENCHMARK"

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
