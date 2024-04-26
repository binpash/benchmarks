#!/bin/bash

unix50() {
    times_file="../outputs/time.res"
    hashed_file="../outputs/hashed.res"
    suite_suffix="bash"
    outputs_dir="../outputs"

    if [ -e "$times_file" ]; then
        echo "skipping unix50"
        return 0
    fi

    if [ ! -d $outputs_dir ]; then
        mkdir "$outputs_dir"
    fi

    scripts_inputs=(
        "1;1.txt"
        "2;1.txt"
        "3;1.txt"
        "4;1.txt"
        "5;2.txt"
        "6;3.txt"
        "7;4.txt"
        "8;4.txt"
        "9;4.txt"
        "10;4.txt"
        "11;4.txt"
        "12;4.txt"
        "13;5.txt"
        "14;6.txt"
        "15;7.txt"
        "16;7.txt"
        "17;7.txt"
        "18;8.txt"
        "19;8.txt"
        "20;8.txt"
        "21;8.txt"
        # "22;8.txt"
        "23;9.1.txt"
        "24;9.2.txt"
        "25;9.3.txt"
        "26;9.4.txt"
        # "27;9.5.txt"
        "28;9.6.txt"
        "29;9.7.txt"
        "30;9.8.txt"
        "31;9.9.txt"
        "32;10.txt"
        "33;10.txt"
        "34;10.txt"
        "35;11.txt"
        "36;11.txt"
    )

    touch "$times_file"
    for script_input in ${scripts_inputs[@]}; do
        IFS=";" read -r -a script_input_parsed <<<"${script_input}"
        script="${script_input_parsed[0]}"
        input="${script_input_parsed[1]}"
        export IN="./input_txt/$input"
        printf -v pad %30s
        padded_script="${script}.sh:${pad}"
        padded_script=${padded_script:0:30}

        outputs_file="${outputs_dir}/${script}.${suite_suffix}.out"
        err_file="${outputs_dir}/${script}.${suite_suffix}.err"

        if [ -e "$outputs_file" ]; then
            echo "skipping $script"
            continue
        fi

        echo "${padded_script}" $({ time ../${script}.sh 2> /dev/null > "$outputs_file"; } 2>&1) | tee -a "$times_file"

        #shasum -a 256 in macOS, /exists --check flag
        sha256sum $outputs_file >>$hashed_file

    done
}

echo "Running unix50"
unix50
