#!/bin/bash

"$MISC_SCRIPT_DIR/grep_and_sleep.sh" 0.5 foo "$test_output_dir/in1" "$test_output_dir/out1"
"$MISC_SCRIPT_DIR/grep_and_sleep.sh" 1 bar "$test_output_dir/in2" "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/grep_and_sleep.sh" 1 bar "$test_output_dir/in3" "$test_output_dir/out3"
"$MISC_SCRIPT_DIR/grep_and_sleep.sh" 1 bar "$test_output_dir/in4" "$test_output_dir/out4"
"$MISC_SCRIPT_DIR/grep_and_sleep.sh" 1 foo "$test_output_dir/out1" "$test_output_dir/out6"