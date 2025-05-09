#!/bin/bash

export file_size=${1?No file size given}

head --bytes="${file_size}" /dev/urandom > test.txt
