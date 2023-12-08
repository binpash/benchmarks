#!/bin/bash
convert_to_seconds() {
    delimiter_min="m"
    delimiter_sec="s"
    delimiter_dec=","

    local num1=$1

    local minutes="${num1%%$delimiter_min*}"
    local seconds="${num1#*$delimiter_min}"
    local seconds="${seconds%%$delimiter_sec*}"
    local seconds_before="${seconds%%$delimiter_dec*}"
    local seconds_dec="${seconds#*$delimiter_dec}"

    local total_sec=$((10#$minutes * 60000 + 10#$seconds_before * 1000 + 10#$seconds_dec))
    echo "$total_sec"
}

time_check(){
    local file1=$1
    local file2=$2

    local cut1=$(awk '{print $3}' "$file1")
    local cut2=$(awk '{print $3}' "$file2")
    local cutname=$(awk '{print $1}' "$file1")

    paste <(echo "$cutname") <(echo "$cut1") <(echo "$cut2") > "./temp/joint_cut.txt"

    while read -r num1 num2 num3; do

        local first=$(convert_to_seconds "$num2")
        local second=$(convert_to_seconds "$num3")

        local diff=$((first - second))   
        
        if [ "$diff" -lt 0 ];then
            diff=$((-diff))
        fi

        if [[ $diff -gt $((first * 20 / 100)) ]];then #does not work for very small times :(
            echo "$num1" "Deviation in excecution times! Make sure everything runs smoothly."
        fi

    done < "./temp/joint_cut.txt"
}
echo "Verifying correct execution"
echo ""

if [ ! -d ./temp ]; then
    mkdir ./temp
fi

size="." #changes for small/large inputs
og_output_time="../for_validation/${size}/for_validation.time"
my_output_time="../outputs/${size}/time.res"

echo "Time comparison"
time_check "$og_output_time" "$my_output_time" 

echo ""

echo "Values comparison"
og_output_avg="../for_validation/${size}/for_validation.avg"
og_output_min="../for_validation/${size}/for_validation.min"
og_output_max="../for_validation/${size}/for_validation.max"

my_output_avg="../outputs/${size}/average.txt"
my_output_min="../outputs/${size}/min.txt"
my_output_max="../outputs/${size}/max.txt"


diff -q "$og_output_avg" "$my_output_avg"
diff -q "$og_output_min" "$my_output_min"
diff -q "$og_output_max" "$my_output_max"