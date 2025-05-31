#!/bin/bash

# Combine terms to create  n-grams (for n=1,2,3) and then count and sort them
# Usage: ./combine.sh <terms > n-grams

p1=$(mktemp -u p1.XXXXXX)
p2=$(mktemp -u p2.XXXXXX)
p3=$(mktemp -u p3.XXXXXX)

mkfifo "$p1" "$p2" "$p3"

bigram() {
	tee "$1" | tail +2 | paste "$1" - | sort
	rm "$1"
}

trigram() {
	tee "$1" | tail +2 | paste "$1" - | tee "$2" | cut -f 1 | tail +3 | paste "$2" - | sort
	rm "$1" "$2"
}

tee >(sort) >(bigram "$p1") >(trigram "$p2" "$p3") > /dev/null