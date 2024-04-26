#!/bin/bash

# exit when any command fails
#set -e

IN=$PASH_TOP/evaluation/benchmarks/mediaconv/input/
OUT=$PASH_TOP/evaluation/benchmarks/mediaconv/output/
IN_NAME=$PASH_TOP/evaluation/benchmarks/mediaconv/input/100G.txt

if [ "$1" == "--small" ]; then
    LOG_DATA_FILES=6
    WAV_DATA_FILES=20
    # TODO pac-n4 machine is down
    JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/jpg.zip
else
    LOG_DATA_FILES=84
    WAV_DATA_FILES=120
    # TODO pac-n4 machine is down
    JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/jpg.zip
fi

if [ ! -d ${IN}/wav ]; then
    # TODO pac-n4 machine is down
    wget http://pac-n4.csail.mit.edu:81/pash_data/wav.zip
    unzip wav.zip && cd wav/
    for f in *.wav; do
        FILE=$(basename "$f")
        for (( i = 0; i <= $WAV_DATA_FILES; i++)) do
            echo copying to $f$i.wav
            cp $f $f$i.wav
        done
    done
    echo "WAV Generated"
fi

if [ ! -d ${IN}/jpg ]; then
    cd ${IN}
    wget $JPG_DATA_LINK
    unzip jpg.zip
    echo "JPG Generated"
    rm -rf ${IN}/jpg.zip
fi

# TODO What's the purpose of this?
source_var() {
  export IN=
}
