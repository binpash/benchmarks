#!/bin/bash

small_input="1M.txt"

if [[ "$1" == "-c" ]]; then
    if [[ "$2" == "--small" ]]; then
        cd small
        rm -rf $small_input
        cd ..
    fi
#    exit
fi

setup_dataset_oneliners() {
    
#     if [ ! -f ./dict.txt ]; then
#         curl -sf 'http://ndr.md/data/dummy/dict.txt' | sort > dict.txt
#         if [ $? -ne 0 ]; then
#             echo 'cannot find dict.txt -- please contact the developers of pash'
#             exit 1
#       fi
# #      "$PASH_TOP/scripts/append_nl_if_not.sh" ./dict.txt
#     fi


    if [ ! -d ./temp ]; then
        mkdir temp
        echo "Downloading original text"
        wget -O ./temp/original.txt https://www.gutenberg.org/cache/epub/59306/pg59306.txt
    fi
#    if [[ "$1" == "--small" ]]; then
#        if [ ! -d ./small ]; then
#            mkdir small
#        fi
#        if [ ! -f ./small/1M.txt ]; then
    touch ./small/1M.txt
    echo "Generating small-size inputs"
    while [ $(wc -c <"./small/1M.txt") -lt 1000000 ] 
    do
        cat ./temp/original.txt >> ./small/1M.txt 
    done               
#        fi
#    fi
    if [[ "$1" == "--large" ]]; then
        if [ ! -d ./large ]; then
            mkdir large
        fi
        if [ ! -f ./large/1G.txt ]; then
            touch ./large/1G.txt
            echo "Generating large-size inputs"
            while [ $(wc -c <"./large/1G.txt") -lt 1000000 ] #choose end goal 
            do
                cat ./temp/original.txt >> ./large/1G.txt 
            done               
        fi
    fi

}

echo "Setting up inputs"
setup_dataset_oneliners