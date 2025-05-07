#!/bin/bash

BENCH_TOP=${BENCH_TOP:-$(git rev-parse --show-toplevel)}
RESOURCES_DIR=${RESOURCES_DIR:-$BENCH_TOP/web-index}/inputs
URL="https://atlas.cs.brown.edu/data"

mkdir -p $RESOURCES_DIR

cp $BENCH_TOP/web-index/stopwords.txt $RESOURCES_DIR

is_small=false
is_min=false
for arg in "$@"; do
    if [ "$arg" = "--small" ]; then
        is_small=true
		suffix="_small"
        break
    fi
	if [ "$arg" = "--min" ]; then
		is_min=true
		suffix="_min"
		break
	fi
done

if [[ ! -d "$RESOURCES_DIR/articles$suffix" ]]; then
	if [[ ! -f "$RESOURCES_DIR/wikipedia$suffix.tar.gz" ]]; then
		if $is_min; then
			echo "Downloading the min dataset."
			# previously was the small dataset
			wget -O $RESOURCES_DIR/wikipedia$suffix.tar.gz "${URL}/wikipedia/input_small/articles.tar.gz" --no-check-certificate
			wget -O $RESOURCES_DIR/index$suffix.txt "${URL}/wikipedia/input_small/index.txt" --no-check-certificate
			echo "Extracting the min dataset."
			tar -xf $RESOURCES_DIR/wikipedia$suffix.tar.gz -C $RESOURCES_DIR
			mv $RESOURCES_DIR/articles $RESOURCES_DIR/articles$suffix
		elif $is_small; then
			# 1gb entries
			echo "Downloading the small dataset."
			wget --no-check-certificate -O $RESOURCES_DIR/wikipedia$suffix.tar.gz "${URL}/wikipedia/wikipedia1g.tar.gz"
			wget --no-check-certificate -O $RESOURCES_DIR/index$suffix.txt "${URL}/wikipedia/index1g.txt"
			echo "Extracting the small dataset."
			tar -xf $RESOURCES_DIR/wikipedia$suffix.tar.gz -C $RESOURCES_DIR
			mv $RESOURCES_DIR/articles1g $RESOURCES_DIR/articles$suffix
		else
			# full dataset
			echo "Downloading the full dataset."
			wget --no-check-certificate -O $RESOURCES_DIR/wikipedia$suffix.tar.gz "${URL}/wikipedia/wikipedia10g.tar.gz"
			wget --no-check-certificate -O $RESOURCES_DIR/index$suffix.txt "$URL}/wikipedia/index10g.txt"
			echo "Extracting the full dataset."
			tar -xf $RESOURCES_DIR/wikipedia$suffix.tar.gz -C $RESOURCES_DIR
			mv $RESOURCES_DIR/articles10g $RESOURCES_DIR/articles$suffix
		fi
	else
		echo "Extracting dataset."
		tar -xf $RESOURCES_DIR/wikipedia$suffix.tar.gz -C $RESOURCES_DIR
	fi
else
	echo "Dataset already exists."
fi
