#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p tmp
mkdir -p result
mkdir -p inputs

/usr/bin/env python3 -c "from sklearn.datasets import fetch_kddcup99; fetch_kddcup99(data_home=\"inputs\", percent10=False, download_if_missing=True)"