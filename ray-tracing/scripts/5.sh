#!/bin/bash

csv_file=$1

cat "$csv_file" | q -H -d, "SELECT MAX(timestamp), MAX(hop) FROM - GROUP BY pathID LIMIT 5"
