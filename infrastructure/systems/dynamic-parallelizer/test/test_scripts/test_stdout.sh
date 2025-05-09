## Test commands that write to stdout
cat "$test_output_dir/in1"
"$MISC_SCRIPT_DIR/sleep_and_grep_stdout.sh" 0.4 "foo" "$test_output_dir/in1"
"$MISC_SCRIPT_DIR/sleep_and_grep_stdout.sh" 0.2 "bar" "$test_output_dir/in1"
"$MISC_SCRIPT_DIR/sleep_and_grep_stdout.sh" 0 "baz" "$test_output_dir/in1"
echo "hi"
pwd
