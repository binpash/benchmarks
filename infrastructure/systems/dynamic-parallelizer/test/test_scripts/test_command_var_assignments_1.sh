export foo=$("$MISC_SCRIPT_DIR/sleep_and_return.sh" 0.35 "hello")
export bar=$("$MISC_SCRIPT_DIR/sleep_and_return.sh" 0.3 "world")
export filename=$("$MISC_SCRIPT_DIR/sleep_and_return.sh" 0.25 "out1")
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.2 "$foo" "$test_output_dir/$filename"
export foo=$("$MISC_SCRIPT_DIR/sleep_and_return.sh" 0.15 "HELLO")
"$MISC_SCRIPT_DIR/sleep_and_cat.sh" 0.1 "$test_output_dir/$filename" "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.05 "$foo $bar" "$test_output_dir/out2"
cat "$test_output_dir/out1" "$test_output_dir/out2"
