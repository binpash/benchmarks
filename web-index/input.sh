#!/bin/bash

BENCH_TOP=${BENCH_TOP:-$(git rev-parse --show-toplevel)}
RESOURCES_DIR=${RESOURCES_DIR:-$BENCH_TOP/web-index}

mkdir -p $RESOURCES_DIR

if [ "$1" = "--small" ]; then
	if [[ ! -f "$RESOURCES_DIR/wikipedia-small.tar.gz" ]]; then
		# 1000 entries
		echo "Downloading the small dataset."
		wget -O $RESOURCES_DIR/wikipedia-small.tar.gz https://atlas-group.cs.brown.edu/data/wikipedia/input_small/articles.tar.gz --no-check-certificate
		wget -O $RESOURCES_DIR/index_small.txt https://atlas-group.cs.brown.edu/data/wikipedia/input_small/index.txt --no-check-certificate
	fi
else
	if [[ ! -f "$RESOURCES_DIR/wikipedia.tar.gz" ]]; then
		# full dataset
		echo "Downloading the full dataset. Caution!! Extracted size >200GB"
		wget -O $RESOURCES_DIR/wikipedia.tar.gz https://atlas-group.cs.brown.edu/data/wikipedia/input/articles.tar.gz --no-check-certificate
		wget -O $RESOURCES_DIR/index.txt https://atlas-group.cs.brown.edu/data/wikipedia/input/index.txt --no-check-certificate
	fi
fi

if [[ ! -d "$RESOURCES_DIR/articles" ]]; then
	if [ "$1" = "--small" ]; then
		# 1000 entries
		echo "Extracting the small dataset."
		tar -xf $RESOURCES_DIR/wikipedia-small.tar.gz -C $RESOURCES_DIR
	else
		# full dataset
		echo "Extracting the full dataset. Caution!! Extracted size >200GB"
		tar -xf $RESOURCES_DIR/wikipedia.tar.gz -C $RESOURCES_DIR
	fi
else
	echo "Did not extract data because of existing data."
	echo "Please rm -r $RESOURCES_DIR/articles manually and rerun this script."
fi

echo "Data is ready."
