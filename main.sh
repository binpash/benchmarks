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
    REPO_TOP="/benchmarks"


    cd "$(dirname "$0")/$BENCHMARK" || exit 1

    # Download dependencies
    ./deps.sh $@ || error "Failed to download dependencies for $BENCHMARK"

    # Fetch inputs
    ./input.sh $@ || error "Failed to fetch inputs for $BENCHMARK"


    if $measure_resources; then
        sudo apt-get install -y autoconf automake libtool build-essential cloc
        cd "$REPO_TOP" || exit 1
        pip install --break-system-packages -r "$REPO_TOP/infrastructure/requirements.txt"
        python3 "$REPO_TOP/infrastructure/run_dynamic.py" $BENCHMARK $@  || error "Failed to run $BENCHMARK"

        rm "$REPO_TOP/infrastructure/target/dynamic_analysis.jsonl"

        cd "$REPO_TOP/infrastructure" || exit 1
        make target/dynamic_analysis.jsonl
        python3 viz/dynamic.py "$REPO_TOP/$BENCHMARK" --text
        cat "$REPO_TOP/$BENCHMARK/benchmark_stats.txt"
        cd "$REPO_TOP/$BENCHMARK" || exit 1
    elif $measure_time; then
        /usr/bin/time -f "Runtime: %E (CPU: %P)" ./run.sh $@ > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
    else
        ./run.sh $@ > "$BENCHMARK.out" 2> "$BENCHMARK.err" || error "Failed to run $BENCHMARK"
    fi

    # Verify output
    ./verify.sh $@ > "$BENCHMARK.hash" || error "Failed to verify output for $BENCHMARK"

    # Cleanup outputs
    ./cleanup.sh $@ 

    if correct; then
        echo "$BENCHMARK [pass]"
    else
        error "$BENCHMARK [fail]"
    fi

    cd - || exit 1
}

main "$@"
