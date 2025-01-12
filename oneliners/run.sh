#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

export BENCHMARK_CATEGORY="oneliners"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

if [[ "$@" == *"--small"* ]]; then
    scripts_inputs=(
        "nfa-regex;1M"
        "sort;1M"
        "top-n;1M"
        "wf;1M"
        "spell;1M"
        "diff;1M"
        "bi-grams;1M"
        "set-diff;1M"
        "sort-sort;1M"
        "uniq-ips;logs-popcount-org"
    )
else
    scripts_inputs=(
        "nfa-regex;1G"
        "sort;3G"
        "top-n;3G"
        "wf;3G"
        "spell;3G"
        "diff;3G"
        "bi-grams;3G"
        "set-diff;3G"
        "sort-sort;3G"
        "uniq-ips;logs-popcount-org"
    )
fi

mkdir -p "outputs"

echo executing oneliners $(date)

for script_input in ${scripts_inputs[@]}
do
    IFS=";" read -r -a parsed <<< "${script_input}"
    script_file="./scripts/${parsed[0]}.sh"
    input_file="./inputs/${parsed[1]}.txt"
    output_file="./outputs/${parsed[0]}.out"

    echo "$script_file"
    export BENCHMARK_INPUT_FILE="$(realpath "$input_file")"
    export BENCHMARK_SCRIPT="$(realpath "$script_file")"
    $BENCHMARK_SHELL "$script_file" "$input_file" > "$output_file"
    echo "$?"
done
