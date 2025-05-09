#!/bin/bash

# Extend the PATH and define key directories
export PATH=$PATH:$HOME/.local/bin
export PASH_SPEC_TOP="${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}"
export PASH_TOP="${PASH_TOP:-$PASH_SPEC_TOP/deps/pash}"

# Define directory for outputs
out_dir="$PASH_SPEC_TOP/report/output/unix_50/full/"
new_out_dir="$PASH_SPEC_TOP/report/output/unix_50/"
script_dir="$PASH_SPEC_TOP/report/benchmarks/unix_50/full"

# Ensure output directory exists
mkdir -p "$out_dir"

# Cleanup specific subdirectories if they exist
if [ -d "$out_dir" ]; then
    rm -rf "$new_out_dir/base"
fi


# Run baseline benchmark with sh-only
echo "Running command: $script_dir/run --target sh-only"
"$script_dir/run" --target sh-only

# Rename output directory to 'base' for the baseline run
mkdir -p "$new_out_dir/base/"
mv "$out_dir/"* "$new_out_dir/base/"

# Run benchmarks with hs-only for each window size
for window in 0 10 20 30 40; do
    # Cleanup specific subdirectories before each run
    if [ -d "$out_dir" ]; then
        rm -rf "$new_out_dir/$window"
    fi

    # Run the benchmark with hs-only target and the specified window
    "$script_dir/run" --target hs-only --window "$window"

    mkdir -p "$new_out_dir/$window/"
    mv "$out_dir/"* "$new_out_dir/$window/"

done
