
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

    echo "$minutes" "$seconds_before" "$seconds_dec"

    local total_sec=$((minutes * 60000 + seconds_before * 1000 + seconds_dec))
    echo "$total_sec"
}

cut1=$(awk '{print $3}' "$file1")
cut2=$(awk '{print $3}' "$file2")

paste <(echo "$cut1") <(echo "$cut2") > "./temp/joint_cut.txt"

while read -r num1_og num2_og; do

    num1=$(convert_to_seconds "$num1_og")
    num2=$(convert_to_seconds "$num2_og")

    #TODO, COMPARE WITH SOME CRITERION? 

    if [[ "$num1" > "$num2" ]]; then
        echo "$num1_og is greater than $num2_og"
    elif [[ "$num1" < "$num2" ]]; then
        echo "$num1_og is less than $num2_og"
    else
        echo "$num1_og is equal to $num2_og"
    fi

done < "./temp/joint_cut.txt"