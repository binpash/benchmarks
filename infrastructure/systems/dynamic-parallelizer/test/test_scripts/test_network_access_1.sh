"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 0.7 "foo" "$test_output_dir/in1" "$test_output_dir/out1"
"$MISC_SCRIPT_DIR/sleep_and_grep.sh" 0.4 "foo" "$test_output_dir/out1" "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/sleep_and_curl.sh" 0.1
