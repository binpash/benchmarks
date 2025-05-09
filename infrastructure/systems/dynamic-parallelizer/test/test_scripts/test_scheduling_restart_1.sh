#!/bin/bash

"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 1 foo "$test_output_dir/in1" "$test_output_dir/out1"
"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 3 bar "$test_output_dir/in2" "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 1 foo "$test_output_dir/out1" "$test_output_dir/out3"
