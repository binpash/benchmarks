# red color
RED='\033[0;31m'
# reset the color
NC='\033[0m'

IN=${BIO4:-$PASH_TOP/benchmarks/bio}
IN_NAME=${IN_N:-input.txt}
if [[ $1 == "-c" ]]; then
    rm -rf *.bam
    rm -rf *.sam
    rm -rf ../output
    exit
fi

PW=${PASH_TOP}/benchmarks/bio/input
echo $PW
mkdir -p $PW
mkdir -p ${PASH_TOP}/benchmarks/bio/output

cat ${IN}/${IN_NAME} |while read s_line;
	do
    
    sample=$(echo $s_line |cut -d " " -f 2);
    if [[ ! -f $sample ]]; then
        pop=$(echo $s_line |cut -f 1 -d " ");
        link=$(echo $s_line |cut -f 3 -d " ");
        wget -O "$PW/$sample".bam  "$link"; ##this part can be adjusted maybe
        # stop after one download
        exit 0
    fi
done;
