#!/bin/bash

BENCH_TOP=${BENCH_TOP:-$(git rev-parse --show-toplevel)}
RESOURCES_DIR=${RESOURCES_DIR:-$BENCH_TOP/report/resources/web-index/}

mkdir -p $RESOURCES_DIR

# echo "Downloading the full dataset. Caution!! Extracted size >200GB"
# wget --no-check-certificate -O $RESOURCES_DIR/wikipedia.tar.gz https://atlas.cs.brown.edu/data/wikipedia/input/articles.tar.gz
# wget --no-check-certificate -O $RESOURCES_DIR/index.txt https://atlas.cs.brown.edu/data/wikipedia/input/index.txt
# tar -xf $RESOURCES_DIR/wikipedia.tar.gz -C $RESOURCES_DIR

echo "Downloading the 100MB dataset"
wget --no-check-certificate -O $RESOURCES_DIR/wikipedia100m.tar.gz https://atlas.cs.brown.edu/data/wikipedia/wikipedia100m.tar.gz
wget --no-check-certificate -O $RESOURCES_DIR/index100m.txt https://atlas.cs.brown.edu/data/wikipedia/index100m.txt
tar -xf $RESOURCES_DIR/wikipedia100m.tar.gz -C $RESOURCES_DIR

echo "Downloading the 1G dataset"
wget --no-check-certificate -O $RESOURCES_DIR/wikipedia1g.tar.gz https://atlas.cs.brown.edu/data/wikipedia/wikipedia1g.tar.gz
wget --no-check-certificate -O $RESOURCES_DIR/index1g.txt https://atlas.cs.brown.edu/data/wikipedia/index1g.txt
tar -xf $RESOURCES_DIR/wikipedia1g.tar.gz -C $RESOURCES_DIR

# echo "Downloading the 10G dataset"
# wget --no-check-certificate -O $RESOURCES_DIR/wikipedia10g.tar.gz https://atlas.cs.brown.edu/data/wikipedia/wikipedia10g.tar.gz
# wget --no-check-certificate -O $RESOURCES_DIR/index10g.txt https://atlas.cs.brown.edu/data/wikipedia/index10g.txt
# tar -xf $RESOURCES_DIR/wikipedia10g.tar.gz -C $RESOURCES_DIR
