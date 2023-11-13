#!/bin/bash

echo "Verifying"

time_check(){
    otherline=$(cat $2)
    counter=0

    for line in $(cat $1)
    do
        echo
    done
}


compare_outputs(){
    if [ "$1" == "--small" ]; then
        size="." #change
    else
        size = "large"
    fi
    og_output_time="../for_validation/${size}/for_validation.time"
    my_output_time="../outputs/${size}/asimpletime.res"

    echo "Time comparison"
    time_check "$og_output_time" "$my_output_time"  #todo: actually compare time. Display warning over some time difference

    echo "Hash comparison"
    og_output_hash="../for_validation/${size}/for_validation.hash"
    my_output_hash="../outputs/${size}/hashed.res"
    diff -q "$og_output_hash" "$my_output_hash"
}

compare_outputs $1