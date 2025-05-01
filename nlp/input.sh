#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1
mkdir -p inputs
cd inputs || exit 1

if [ ! -f ./book_links.txt ]; then
    wget -O book_links.txt "https://atlas-group.cs.brown.edu/data/gutenberg/books.txt"
    if [ ! -f book_links.txt ]; then
        echo "Failed to download book_links.txt"
        exit 1
    fi
fi

if [ ! -f ./genesis ]; then
    curl -sf https://atlas-group.cs.brown.edu/data/gutenberg/8/0/0/8001/8001.txt > genesis
fi 

if [ ! -f ./exodus ]; then
    curl -sf https://atlas-group.cs.brown.edu/data/gutenberg/3/3/4/2/33420/33420-0.txt > exodus
fi

for arg in "$@"; do
    if [[ "$arg" == "--small" ]]; then
        if [ ! -e ./pg-small ]; then
            mkdir pg-small
            cd pg-small || exit 1
            book_count=3000

            head -n $book_count ../book_links.txt | while IFS= read -r line
            do
                full_url="https://atlas-group.cs.brown.edu/data/gutenberg/${line}"
                echo "Downloading $full_url"
                wget -q "$full_url"
            done
            cd ..
            exit 0
        fi
    elif [[ "$arg" == "--min" ]]; then
        if [ ! -e ./pg-min ]; then
            mkdir pg-min
            cd pg-min || exit 1
            book_count=1

            head -n $book_count ../book_links.txt | while IFS= read -r line
            do
                full_url="https://atlas-group.cs.brown.edu/data/gutenberg/${line}"
                echo "Downloading $full_url"
                wget -q "$full_url"
            done

            cd ..
            exit 0
        fi
    fi
done

if [ ! -e ./pg ]; then
    mkdir pg
    cd pg || exit 1

    while IFS= read -r line
    do
        full_url="https://atlas-group.cs.brown.edu/data/gutenberg/${line}"
        echo "Downloading $full_url"
        wget -q "$full_url"
    done < ../book_links.txt

    cd ..
fi
