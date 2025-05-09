#!/bin/sh

for test in $(cat script_names.txt); do
    sh $test
done
