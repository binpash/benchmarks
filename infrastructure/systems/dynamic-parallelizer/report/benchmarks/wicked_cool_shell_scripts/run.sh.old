#!/bin/sh
#set -x

export SYSADMINDIR="${PASH_SPEC_TOP}/report/benchmarks/wicked_cool_shell_scripts"
export OUTBASE="${PASH_SPEC_TOP}/report/output/wicked_cool_shell_scripts"

if [ -z window ]
then
    echo please set window var
    exit 1
fi

pashcmd="/srv/hs/pash-spec.sh --window $window"

if [ -z logging ]
then
    pashcmd="${pashcmd} -d2"
fi

run() {
    local job=$1
    cd $SYSADMINDIR/$job

    echo running sh
    /usr/bin/time -f '%e' -o result_sh_time sh run_sh.sh > result_sh_log 2>&1
    echo running hs
    /usr/bin/time -f '%e' -o result_hs_time $pashcmd run_hs.sh > result_hs_log 2>&1
    echo verifying
    sh verify.sh > result_error
    cat result_error
    mkdir -p $OUTBASE/$job
    mv result_* $OUTBASE/$job
}

echo --- running 27 ---
run 27
echo --- running 28 ---
run 28
echo --- running 29 ---
run 29
# TODO 35 currently hangs
#echo --- running 35 ---
#run 35
# TODO 37 df output will be different
#echo --- running 37 ---
#run 37
# TODO 40 verify fails due to shadow file owned by another gid
#echo --- running 40 ---
#run 40
# TODO 41 as it is an interactive focused script
#echo --- running 41 ---
#run 41
# TODO 42 init fails due to shadow file owned by another gid
#echo --- running 42 ---
#run 42
echo --- running 43 ---
run 43
echo --- running 48 ---
run 48
echo --- running 50 ---
run 50
# TODO 51 proposed removal due to timestamps involved
#echo --- running 51 ---
#run 51
echo --- running 52 ---
run 52
echo --- running 102 ---
run 102

echo all done
