#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import json
from subprocess import check_output
from collections import Counter
import sys

from all_scripts import get_all_scripts
from syntax_analysis import parse_shell_script, count_nodes
from project_root import get_project_root

run(['bash', *sys.argv[1:]])
