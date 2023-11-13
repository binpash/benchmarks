
#!/bin/bash    
size='.'
file1="../for_validation/${size}/for_validation.time"
file2="../outputs/${size}/asimpletime.res"

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

    local total_sec=$((minutes * 60000 + seconds_before * 1000 + seconds_dec))
    echo "$total_sec"
}

cut1=$(awk '{print $3}' "$file1")
cut2=$(awk '{print $3}' "$file2")
cutname=$(awk '{print $1}' "$file1")

paste <(echo "$cutname") <(echo "$cut1") <(echo "$cut2") > "./temp/joint_cut.txt"

while read -r num1 num2 num3; do

    first=$(convert_to_seconds "$num2")
    second=$(convert_to_seconds "$num3")

    diff=$((first - second))
    if [ "$diff" -lt 0 ];then
        diff=$((-diff))
    fi

    if [[ $diff -gt $((first * 2 / 100)) ]];then #does not work for very small times :(
        echo "$num1" "Deviation in excecution times! Make sure everything runs smoothly."
    fi

done < "./temp/joint_cut.txt"