#!/bin/bash

sleep $1
grep $2 $3 > $4
echo "grepping $1 $2"
