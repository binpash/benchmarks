#!/bin/bash

csv_file=$1

cat "$csv_file" | q -H -d, "SELECT * FROM - WHERE pathID = 20613314"