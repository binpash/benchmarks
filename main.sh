error()
{
    echo "Error: $1" > /dev/stderr
    exit 1
}

correct() { [ "$(cat $BENCHMARK.hash | cut -d' ' -f 2 | grep -c 1)" -eq 0 ]; }

main()
{
    cd "$(dirname "$0")/$BENCHMARK" || exit 1
    # Download dependencies
    ./deps.sh $@ || error "Failed to download dependencies for $BENCHMARK"

    # Fetch inputs
    ./inputs.sh $@  || error "Failed to fetch inputs for $BENCHMARK"

    # Run benchmark
    ( ./run.sh $@ > $BENCHMARK.out 2> $BENCHMARK.err ) || error "Failed to run $BENCHMARK"

    # Verify output
    ./verify.sh $@ > $BENCHMARK.hash || error "Failed to verify output for $BENCHMARK"

    if correct; then
        echo "$BENCHMARK [pass]"
    else
        error "$BENCHMARK [fail]"
    fi

    cd - || exit 1
}

export BENCHMARK="$1"
shift

main $@
