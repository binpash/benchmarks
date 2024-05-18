#!/bin/bash

# exit when any command fails
#set -e

# Download inputs
# 						~ Use two sources: University (fast) 
# 					             			Long term storage (slow)
# 						~ And two sizes:  Small & quick
# 	             Full size

IN=inputs
OUT=outputs
# IN=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/
# OUT=$PASH_TOP/evaluation/benchmarks/dependency_untangling/output/
# IN_NAME=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/100G.txt

if [ "$1" == "-c" ]; then
    rm -rf ${IN}/pcap_data
    rm -rf ${IN}/pcaps
    rm -rf ${OUT}
    exit 
fi

setup_dataset() {
  if [ "$1" == "--small" ]; then
      PCAP_DATA_FILES=1
  else
      PCAP_DATA_FILES=15
  fi                                                                       
  
  # download the initial pcaps to populate the whole dataset
  if [ ! -d ${IN}/pcap_data ]; then
    # wget https://ita.ee.lbl.gov/traces/clarknet_access_log_Aug28.gz
    # wget https://ita.ee.lbl.gov/traces/clarknet_access_log_Sep4.gz

    #   wget http://pac-n4.csail.mit.edu:81/pash_data/pcaps.zip
    #   unzip pcaps.zip
    #   rm pcaps.zip
    mkdir ${IN}/pcap_data/
    # Make it 20G
    for item in ${IN}/pcaps/*;
    do
        name=$(basename $item)
        filename_final="${IN}/pcap_data/pcap_${name}"
        for ((i = 0; i < 20; i++)); do
            cat $item >> $filename_final
        done
        # Convert the compresed ASCII file into pcap format
        tcpdump -r "$filename_final" -w "$filename_final.pcap"
    done
    #   for ((file_num = 0; file_num < 10; file_num++)); do
    #     filename="${IN}/pcaps/alias_${file_num}.sh"
    #     filename_final="${IN}/pcap_data/pcap_${file_num}.sh"
    #     touch $filename_final
    #     echo $filename_final
    #     for ((i = 0; i < 100; i++)); do
    #         cat $filename >> $filename_final
    #     done
    #   done
    echo "Pcaps Generated"
  fi 
}

source_var() {
  export IN=
}

setup_dataset