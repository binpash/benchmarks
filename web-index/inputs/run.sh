source_var() {
    export WEB_INDEX_DIR=../data
    export WIKI=../data/
}


web-index(){
    times_file="seq.res"
    outputs_suffix="seq.out"
    outputs_dir="../outputs"
    if [ -e "web-index/${times_file}" ]; then
        echo "skipping web-index/${times_file}"
        return 0
    fi

    if [[ ! -d "$outputs_dir" ]]; then
        mkdir -p "$outputs_dir"
    fi

    touch "$times_file"
    outputs_file="${outputs_dir}/web-index.${outputs_suffix}"
    echo web-index.sh: $({ time ../web-index.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"
}

