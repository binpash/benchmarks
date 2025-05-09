#!/bin/bash

export PATH=$PATH:$HOME/.local/bin
export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
export BIODIR="${PASH_SPEC_TOP}/report/benchmarks/bio4-tryoverhead"
export OUTBASE="${PASH_SPEC_TOP}/report/output/bio4-$SIZE"
export INPUT_LIST=$SIZE/list
export IN="${PASH_SPEC_TOP}/report/resources/bio4/"

cd "$BIODIR"

mkdir -p "$OUTBASE"

hs_r() {
	local size=$1
	local window=$2
	local log=$3

	echo Running hS test for $size with window $window
	/usr/bin/time -f '%e' -o "$OUTBASE/try_time" sh fxbio4-try.sh $size &> "$OUTBASE/try_log"
	md5sum $BIODIR/output/* > "$OUTBASE/try_hash"

	rm -rf output
	diff "$OUTBASE/try_hash" "$OUTBASE/sh_hash" > "$OUTBASE/error"
}

sh_r() {
	local size=$1

	echo Running sh baseline for $size
	/usr/bin/time -f '%e' -o "$OUTBASE/sh_time" sh fxbio4.sh $size &> "$OUTBASE/sh_log"
	md5sum $BIODIR/output/* > "$OUTBASE/sh_hash"

	rm -rf output
}

usage() {
	echo "Usage: $0 [--window WINDOW_SIZE] [--target TARGET] [--log LOG_OPTION]"
	echo "  --window WINDOW_SIZE  Window size to run hs with (default: 5)"
	echo "  --target TARGET       Target to run: hs-only, sh-only, or both"
	echo "  --log LOG_OPTION      Whether to enable logging for hs: enable or disable (default: enable)"
	exit 1
}

window=5
target=""
log="enable"

while [ $# -gt 0 ]; do
	case "$1" in
		--window)
			window="$2"
			shift 2
			;;
		--target)
			target="$2"
			shift 2
			;;
		--log)
			log="$2"
			shift 2
			;;
		*)
			usage
			;;
	esac
done

if [ -z "$target" ]; then
	echo "Error: --target argument is required"
	usage
fi

run_hs=false
run_sh=false

case "$target" in
    hs-only)
        run_hs=true
        ;;
    sh-only)
        run_sh=true
        ;;
    both)
        run_hs=true
        run_sh=true
        ;;
esac

if $run_sh
then
sh_r $INPUT_LIST
fi

if $run_hs
then
hs_r $INPUT_LIST $window $log
fi
