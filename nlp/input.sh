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

BOOKS_URL="https://atlas-group.cs.brown.edu/data/gutenberg/books.txt"
BASE_URL="https://atlas-group.cs.brown.edu/data/gutenberg/"

# Download the list of book urls
if [ ! -f ./book_links.txt ]; then
    wget -q -O book_links.txt "$BOOKS_URL"
    if [ ! -f book_links.txt ]; then
        echo "Failed to download book_links.txt"
        exit 1
    fi
fi

if [ ! -f ./genesis ]; then
    # original link: https://www.gutenberg.org/cache/epub/8001/pg8001.txt
    curl -sf https://atlas-group.cs.brown.edu/data/gutenberg/8/0/0/8001/8001.txt > genesis
    "$PASH_TOP/scripts/append_nl_if_not.sh" genesis
fi 

if [ ! -f ./exodus ]; then
    # original link: https://www.gutenberg.org/files/33420/33420-0.txt
    curl -sf https://atlas-group.cs.brown.edu/data/gutenberg/3/3/4/2/33420/33420-0.txt > exodus
    "$PASH_TOP/scripts/append_nl_if_not.sh" exodus
fi

if [ ! -e ./pg ]; then
    mkdir pg
    cd pg
    book_count=10
    if [[ "$1" == "--full" ]]; then
        book_count=1000
    fi

    head -n $book_count ../book_links.txt | while IFS= read -r line
    do
        full_url="${BASE_URL}${line}"
        echo "Downloading $full_url"
        wget -q "$full_url"
    done

    for f in *.txt; do
        append_nl_if_not $f
    done
    cd ..
fi
