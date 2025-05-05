#!/bin/bash

# shell script to run verify.py
parsed_args=()
for arg in "$@"; do
    case "$arg" in
        --small|--min)
            parsed_args+=("$arg")
            ;;
    esac
done
# run the Python script
python3 validate.py "${parsed_args[@]}"

# check if the script ran successfully
echo "sklearn $?"