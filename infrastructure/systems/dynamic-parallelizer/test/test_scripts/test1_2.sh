"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 0.05 "foo" "$test_output_dir/in1" "$test_output_dir/out1"
"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 0.15 "foo" "$test_output_dir/out1" "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 0.25 "foo" "$test_output_dir/out2" "$test_output_dir/out3"
pwd
