#!/bin/bash

if [[ "$1" == "-c" ]]; then
    rm -rf ./input_txt
    exit
fi

setup_dataset_oneliners() {

    if [ ! -d ./temp ]; then
        mkdir temp
        echo "Downloading original text"
        wget -O ./temp/original.txt https://www.gutenberg.org/cache/epub/59306/pg59306.txt
    fi

    if [ ! -d ./input_txt ]; then
        mkdir input_txt
    fi

    echo "Generating artificial inputs"

    touch ./input_txt/1M.txt ./input_txt/10M.txt ./input_txt/100M.txt ./input_txt/1G.txt 
    
    while [ $(wc -c <"./input_txt/1M.txt") -lt 1000000 ] 
    do
        cat ./temp/original.txt >> ./input_txt/1M.txt 
    done               

    while [ $(wc -c <"./input_txt/10M.txt") -lt 10000000 ] 
    do
        cat ./input_txt/1M.txt >> ./input_txt/10M.txt 
    done    

    while [ $(wc -c <"./input_txt/100M.txt") -lt 100000000 ] 
    do
        cat ./input_txt/10M.txt >> ./input_txt/100M.txt 
    done    

    while [ $(wc -c <"./input_txt/1G.txt") -lt 1000000000 ] 
    do
        cat ./input_txt/100M.txt >> ./input_txt/1G.txt 
    done    

    if [ ! -f ./input_txt/dict.txt ]; then
        curl -sf 'https://atlas-group.cs.brown.edu/data/dummy/dict.txt' | sort > ./input_txt/dict.txt
        if [ $? -ne 0 ]; then
            echo 'cannot find dict.txt'
            exit 1
        fi
    fi

}

echo "Setting up inputs"
setup_dataset_oneliners 