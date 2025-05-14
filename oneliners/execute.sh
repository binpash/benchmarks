#!/bin/bash

SUITE_DIR="$(realpath "$(dirname "$0")")"
export SUITE_DIR
cd "$SUITE_DIR" || exit 1

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/oneliners"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/inputs"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"
export TIMEFORMAT=%R

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_CATEGORY="oneliners"
size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done
export LC_ALL=C
if [[ $size == "small" ]]; then
    scripts_inputs=(
        "nfa-regex;10M"
        "sort;30M"
        "top-n;30M"
        "wf;30M"
        "spell;30M"
        "diff;30M"
        "bi-grams;30M"
        "set-diff;30M"
        "sort-sort;30M"
        "uniq-ips;logs-popcount-org"
        "log-search;log-search"
        "opt-parallel;chessdata_small"

    )
    chess_input="$input_dir/chessdata_min"
elif [[ $size == "min" ]]; then
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
        "log-search;log-search"
        "opt-parallel;chessdata_min"
    )
    chess_input="$input_dir/chessdata_small"
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
        "log-search;log-search"
        "opt-parallel;chessdata"
    )
    chess_input="$input_dir/chessdata"
fi

export LC_ALL=C

echo "executing oneliners $(date)"

for script_input in "${scripts_inputs[@]}"
do
    if [[ "$script_input" != "opt-parallel"* ]]; then
        IFS=";" read -r -a parsed <<< "${script_input}"
        script_file="$scripts_dir/${parsed[0]}.sh"
        input_file="$input_dir/${parsed[1]}.txt"
        output_file="$outputs_dir/${parsed[0]}.out"

        echo "$script_file"
        BENCHMARK_INPUT_FILE="$(realpath "$input_file")"
        export BENCHMARK_INPUT_FILE

        BENCHMARK_SCRIPT="$(realpath "$script_file")"
        export BENCHMARK_SCRIPT
        
        $KOALA_SHELL "$script_file" "$input_file" > "$output_file"
        echo "$?"
    else
        IFS=";" read -r -a parsed <<< "${script_input}"
        script_file="$scripts_dir/${parsed[0]}.sh"
        echo "$script_file"
        export BENCHMARK_INPUT_FILE="${chess_input}"
        BENCHMARK_SCRIPT="$(realpath "$script_file")"
        export BENCHMARK_SCRIPT
        $KOALA_SHELL "$script_file" "" > "${outputs_dir}/opt-parallel_$size.out"
        echo "$?"
    fi
done

