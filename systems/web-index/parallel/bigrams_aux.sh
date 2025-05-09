#!/bin/bash

( mkfifo s2 > /dev/null ) ;
( mkfifo s3 > /dev/null ) ;

sed '$d' s2 > s3 &
tee s2 |
    tail +2 |
    paste s3 -
rm s2
rm s3
