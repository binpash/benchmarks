# Tests that backward dependencies get scheduled correctly

"$MISC_SCRIPT_DIR/sleep_and_cat.sh" 0.2 "$test_output_dir/in1" "$test_output_dir/out1"
"$MISC_SCRIPT_DIR/sleep_and_cat.sh" 0.1 "$test_output_dir/out1" > "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/sleep_and_cat.sh" 0 "$test_output_dir/in2" "$test_output_dir/in1"
