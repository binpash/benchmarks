exec 3> "$test_output_dir/out1"
exec 4> "$test_output_dir/out2"
echo output1 >&3
echo output2 >&4

