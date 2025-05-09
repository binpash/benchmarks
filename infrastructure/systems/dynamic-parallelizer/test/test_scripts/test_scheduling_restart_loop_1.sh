#!/bin/bash

for i in $(seq 1 3); do
    "$MISC_SCRIPT_DIR/grep_and_sleep.sh" 1 foo "$test_output_dir/in$i" "$test_output_dir/temp$i"
    "$MISC_SCRIPT_DIR/grep_and_sleep.sh" 3 bar "$test_output_dir/in_2_$i" "$test_output_dir/out_2_$i"
    "$MISC_SCRIPT_DIR/grep_and_sleep.sh" 1 foo "$test_output_dir/temp$i" "$test_output_dir/out$i"
done
