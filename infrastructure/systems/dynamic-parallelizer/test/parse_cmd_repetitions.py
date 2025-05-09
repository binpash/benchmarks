#!/bin/env python3
import sys
import re

# 1: Repetitions to check --single:   total repetitions
#                         --detailed: repetitions dfor each command
# 2: Log file produced by orch. 

def check_detailed_repetitions(lines):
    REGEX = re.compile(r"DEBUG\|.*Executions\|")
    lines = list(filter(REGEX.match, lines))
    lines = [line.split("|")[3].split(",")[1] for line in lines]
    print(" ".join(lines))
    
def check_total_repetitions(lines):
    REGEX = re.compile(r"DEBUG\|.*TotalExec\|")
    lines = list(filter(REGEX.match, lines))
    lines = lines[0].split("|")[3]
    print(lines)

with open(sys.argv[2], "r", encoding="UTF8") as f:
    lines = f.read().split("\n")
    
if sys.argv[1] == "--total":
    check_total_repetitions(lines)
elif sys.argv[1] == "--detailed":
    check_detailed_repetitions(lines)
else:
    assert False