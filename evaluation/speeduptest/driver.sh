#!/bin/bash
# Generates an executable script for running all parallel experiments
# ./driver.sh [parallelism_factor]

set -e

if [[ $(uname) == 'Darwin' ]]; then
  DIFF_STAT="wc -l"
  OUT_DIR="."
  IN=${IN:-../scripts/input/10M.txt}
else
  DIFF_STAT="diffstat"
  OUT_DIR=/dev/shm
  IN=${IN:-../scripts/input/100M.txt}
fi

PRE="dish"
CREATE="touch" # or mkfifo
CPUs=${1:-$(nproc)}
OUT1=${OUT:-./out1.txt}
OUT2=${OUT:-./out2.txt}
SEQ=${2:-"./seq-grep"} # grep has 3 levels

# divide by number of chunks in AWK
echo dividing input to $CPUs chunks
total_size=$(wc -l $IN | awk -F " " '{print $1}') 
chunk_size=$((total_size / CPUs))
split -l $chunk_size $IN $PRE-chunk-

find .  -maxdepth 1 -type p -delete

echo '#!/bin/bash' > $PRE-execute.sh
chmod +x ./$PRE-execute.sh
echo "# This script is auto-generated by driver.sh" >> $PRE-execute.sh
echo "#seq script: time (cat $IN | $SEQ > $OUT)" >> $PRE-execute.sh

# echo "set -x" >> $PRE-execute.sh

echo creating $CPUs channels
counter=0
for chunk in $PRE-chunk-*; do
  # echo "mkfifo $PRE-channel-$((counter++))" >> $PRE-execute.sh
  if [[ $CREATE == 'mkfifo' ]]; then
    $CREATE $OUT_DIR/$PRE-channel-$((counter++))
  fi
done

counter=0
for chunk in $PRE-chunk-*; do
  if [[ $CREATE == 'touch' ]]; then
    # echo 'Channel is persistent file, using `>` to create it'
    echo "cat $chunk | $SEQ > $OUT_DIR/$PRE-channel-$((counter++)) &" >> $PRE-execute.sh
  else
    # echo 'Channel is FIFO, using `>>` to append to it'
    echo "cat $chunk | $SEQ >> $OUT_DIR/$PRE-channel-$((counter++)) &" >> $PRE-execute.sh
  fi
done

# #FIXME: bash doesn't expand `*` in _numberic_ order (1, 10, 2..) affecting cat
# echo cat '$OUT_DIR/$PRE-channel-* >>' $OUT2 >> $PRE-execute.sh 
# # echo 'wait' >> $PRE-execute.sh 

echo 'wait' >> $PRE-execute.sh 

counter=0
args=""
for chunk in $PRE-chunk-*; do
  args="$args $OUT_DIR/$PRE-channel-$((counter++))"
done
echo cat $args '>' $OUT2 >> $PRE-execute.sh 

echo Sequential Timing:
time (cat $IN | $SEQ > $OUT1)

echo Parallel Timing: 
time ./$PRE-execute.sh

echo Result Diff:
diff $OUT1 $OUT2 | $DIFF_STAT

find .  -maxdepth 1 -type p -delete
rm $OUT_DIR/$PRE-channel-*