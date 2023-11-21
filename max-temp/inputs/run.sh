if [[ "$1" == "-c" ]]; then
    rm -r ../outputs/*
    return 0
fi

max-temp(){
    outputs_dir="../outputs"
    times_file="${outputs_dir}/time.res"
    hashed_file="${outputs_dir}/hashed.res"
    outputs_suffix=".out"
    local script=$1

    if [ -e "${times_file}" ]; then
        echo "skipping max-temp"
        return 0
    fi 

    if [ ! -d ../data ]; then
        mkdir ../data
    fi
    
    mkdir -p "$outputs_dir"
    touch "$times_file"
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    echo "${script}.sh:" $({ time ../${script}.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"

}

if [[ $1 == "--pre" ]]; then
    script="analytics-preprocess"
elif [[ $1 == "--pro" ]]; then
    script="analytics-process"
elif [[ $1 == "--full" ]] || [[ $1 == "" ]] ; then
    script="temp-analytics"
else
    echo "Invalid choice! Choices are --pre (for preprocessing), --pro (for processing), --full (for everything)"
    return 0
fi

max-temp $script