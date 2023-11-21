#!/bin/bash
#https://www1.ncdc.noaa.gov/pub/data/noaa/ works!
IN=${IN:-'https://atlas-group.cs.brown.edu/data/noaa/ '}

sed "s;^;$IN;" |
    sed 's;$;/;' |
    xargs -r -n 1 curl -s |
    grep gz |
    tr -s ' \n' |
    cut -d ' ' -f9 |
    sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
    sed  "s;^;$IN;" |
    xargs -n1 curl -s |
    gunzip