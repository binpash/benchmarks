# Tests speculation of completely independent cmds

"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.4 "hello0" "$test_output_dir/out0"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.36 "hello1" "$test_output_dir/out1"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.32 "hello2" "$test_output_dir/out2"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.28 "hello3" "$test_output_dir/out3"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.24 "hello4" "$test_output_dir/out4"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.20 "hello5" "$test_output_dir/out5"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.16 "hello6" "$test_output_dir/out6"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.12 "hello7" "$test_output_dir/out7"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.08 "hello8" "$test_output_dir/out8"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.04 "hello9" "$test_output_dir/out9"
"$MISC_SCRIPT_DIR/sleep_and_echo.sh" 0.0 "hello10" "$test_output_dir/out10"
