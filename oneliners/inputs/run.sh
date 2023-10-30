#!/bin/bash

#have to run with --small (or --large)!

if [[ "$1" == "-c" ]]; then
    rm -r ../outputs/*
    return 0
fi

source_var() {
    if [[ "$1" == "--large" ]]; then
        export IN="./large/$2"
    else
        export IN="./small/$2"
        export dict="./small/dict.txt"
    fi
}

oneliners(){
    smpl_times_file="../outputs/asimpletime.res"
    smpl_hashed_file="../outputs/hashed.res"
    smpl_outputs_suffix="simple.out"
    outputs_dir="../outputs"
        
    if [ -e "$smpl_times_file" ]; then
        echo "skipping oneliners"
        return 0
    fi
    
    if [ ! -d $outputs_dir ]; then
        mkdir "$outputs_dir"
    fi 

  scripts_inputs=(
      "nfa-regex;1M.txt"
      "sort;1M.txt"
      "top-n;1M.txt"
      "wf;1M.txt"
      "spell;1M.txt"
      "diff;1M.txt"
      "bi-grams;1M.txt"
      "set-diff;1M.txt"
      "sort-sort;1M.txt"
#      "shortest-scripts;all_cmdsx100.txt"
 )

    touch "$smpl_times_file"
    echo executing one-liners $(date) | tee -a "$smpl_times_file"
    echo '' >> "$smpl_times_file"

    for script_input in ${scripts_inputs[@]}
    do
        IFS=";" read -r -a script_input_parsed <<< "${script_input}"
        script="${script_input_parsed[0]}"
        input="${script_input_parsed[1]}"
        source_var $1 $input
        printf -v pad %30s                  #pad?
        padded_script="${script}.sh:${pad}"
        padded_script=${padded_script:0:30}

        smpl_outputs_file="${outputs_dir}/${script}.${smpl_outputs_suffix}"

        echo "${padded_script}" $({ time ../${script}.sh > "$smpl_outputs_file"; } 2>&1) | tee -a "$smpl_times_file" 
        
        #shasum -a 256 in macOS, /exists --check flag
        sha256sum $smpl_outputs_file >> $smpl_hashed_file 

    done

}

oneliners $1