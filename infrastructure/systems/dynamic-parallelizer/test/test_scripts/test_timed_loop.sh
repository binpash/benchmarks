for i in 1 2 3 4 5; do
  "$MISC_SCRIPT_DIR/sleep_and_echo.sh" 1 "$i" "$test_output_dir/out${i}"
done
