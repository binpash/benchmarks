#!/bin/bash

append_nl_if_not() {
    if [ -z "$1" ]; then
        echo "No file argument given!"
        exit 1
    else
        if [ ! -f "$1" ]; then
            echo "File $1 doesn't exist!"
            exit 1
        else
            tail -c 1 "$1" | od -ta | grep -q nl
            if [ $? -eq 1 ]
            then
            echo >> "$1"
            fi
        fi
    fi
}


if [ ! -f ./genesis ]; then
    curl -sf https://www.gutenberg.org/cache/epub/8001/pg8001.txt > genesis
    append_nl_if_not genesis
fi 

if [ ! -f ./exodus ]; then
    curl -sf https://www.gutenberg.org/files/33420/33420-0.txt > exodus
    append_nl_if_not exodus
fi

if [ ! -e ./pg ]; then
    mkdir pg
    cd pg
    book_count=10
    if [[ "$1" == "--full" ]]; then
        book_count=1000
    fi

    echo $PWD

    if [ ! -f ../book_txt_links.txt ]; then 
        echo "No link file found!" 
        exit 1
    fi

    head -n $book_count ../book_txt_links.txt | while IFS= read -r line
    do
        echo "Downloading $line"
        # Your code here
        wget -q "$line"
    done

    for f in *.txt; do
        append_nl_if_not $f
    done
    cd ..
fi
