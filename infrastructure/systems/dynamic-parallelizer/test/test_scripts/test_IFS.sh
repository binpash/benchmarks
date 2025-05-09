IFS='
'
echo "line1:line1::line1" > $test_output_dir/data
echo "line2 line2" >> $test_output_dir/data
echo ":line3" >> $test_output_dir/data

for x in $(cat $test_output_dir/data); do
    echo $x
done

unset IFS
for x in $(cat $test_output_dir/data); do
    echo $x
done

IFS=':
'
for x in $(cat $test_output_dir/data); do
    echo $x
done

unset IFS
for x in $(cat $test_output_dir/data); do
    echo $x
done
