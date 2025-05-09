#!/bin/bash
# Match complex regular-expression over input

tr A-Z a-z < $1 | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'