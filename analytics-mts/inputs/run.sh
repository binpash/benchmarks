#!/bin/bash

if [[ "$1" == "-c" ]]; then
    rm -r ../outputs/*
    return 0
fi

source_var() {
  if [[ "$1" == "--small" ]]; then
    export IN="../data/in_small.csv"
  else
    export IN="../data/in.csv"
  fi    
}

analytics-mts(){
    outputs_suffix="out"
    outputs_dir="../outputs"
    times_file="${outputs_dir}/time.res"
    hashed_file="${outputs_dir}/hashed.res"

    if [[ ! -d "$outputs_dir" ]]; then
        mkdir -p "$outputs_dir"
    fi
    
    if [ -e "${times_file}" ]; then
        echo "skipping analytics-mts"
        return 0
    fi

    touch "$times_file"
#    echo executing MTS analytics $(date) | tee -a "$times_file"
    ## FIXME 5.sh is not working yet
    for number in `seq 4`
    do
        script="${number}"
        
        printf -v pad %20s
        padded_script="${script}.sh:${pad}"
        padded_script=${padded_script:0:20}
        # select the respective input
        source_var $1
        outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

        echo "${padded_script}" $({ time ../${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
        #sha256sum $outputs_file >> $hashed_file 

    done
}

source_var $1
analytics-mts