#!/bin/bash

#have to run with --small (or --large)!

if [[ "$1" == "-c" ]]; then
    rm -r ../outputs/*
    return 0
fi

oneliners(){
    smpl_times_file="../outputs/time.res"
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
      "nfa-regex;10M.txt"
      "sort;1G.txt"
      "top-n;100M.txt"
      "wf;1G.txt"
      "spell;100M.txt"
      "diff;100M.txt"
      "bi-grams;100M.txt"
      "set-diff;1G.txt"
      "sort-sort;100M.txt"
  )

    touch "$smpl_times_file"
    for script_input in ${scripts_inputs[@]}
    do
        IFS=";" read -r -a script_input_parsed <<< "${script_input}"
        script="${script_input_parsed[0]}"
        input="${script_input_parsed[1]}"
        export IN="./input_txt/$input"
        printf -v pad %30s                  
        padded_script="${script}.sh:${pad}"
        padded_script=${padded_script:0:30}

        smpl_outputs_file="${outputs_dir}/${script}.${smpl_outputs_suffix}"

        if [ -e "$smpl_outputs_file" ]; then
            echo "skipping $script"
            continue
        fi

        echo "${padded_script}" $({ time ../${script}.sh > "$smpl_outputs_file"; } 2>&1) | tee -a "$smpl_times_file" 
        
        #shasum -a 256 in macOS, /exists --check flag
        sha256sum $smpl_outputs_file >> $smpl_hashed_file 

    done

}
echo "Running onliners"
oneliners