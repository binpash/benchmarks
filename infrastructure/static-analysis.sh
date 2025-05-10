#!/bin/bash

set -e
set -o pipefail

REPO_TOP=$(git rev-parse --show-toplevel)
OUTPUT_PATH="$1"
shift
pushd "$REPO_TOP"/infrastructure

sudo apt-get update && apt-get install -y \
    git \
    python3 python3-pip python3-venv \
    build-essential libtool m4 automake cloc
python3 -m venv venv && . venv/bin/activate && pip install -r requirements.txt


checked_command() {
    local description="$1"
    shift
    echo "Generating $description" >&2
    if ! "$@"; then
	echo "Error generating $description, static analysis failed" >&2
	exit 1
    fi
}

echo "Creating target directory (if it doesn't exist)"
mkdir -p target plots

checked_command "target/scripts_to_benchmark.csv" python3 scripts_to_benchmark.py | sort > target/scripts_to_benchmark.csv

checked_command "target/lines_of_code.csv" python3 count_lines_of_code.py | sort > target/lines_of_code.csv

checked_command "target/nodes_in_scripts.csv" python3 count_nodes_in_scripts.py | sort > target/nodes_in_scripts.csv

echo "Generating target/shellmetrics.sh"
if ! wget --quiet -O target/shellmetrics.sh https://raw.githubusercontent.com/shellspec/shellmetrics/b3bfff2af6880443112cdbf2ea449440b30ab9b0/shellmetrics; then
  echo "Error downloading target/shellmetrics.sh"
  exit 1
fi

checked_command "target/shellmetrics.sh" chmod a+x target/shellmetrics.sh

checked_command "target/cyclomatic.csv" python3 get_cyclomatic.py | sort > target/cyclomatic.csv

checked_command "plots/koala-stx-analysis.pdf" python3 viz/syntax.py plots

mv plots/koala-stx-analysis.pdf "$OUTPUT_PATH"

echo '---------------------------------------------------------------------------------------------------------------------------------------'
echo "Static analysis complete, all targets generated successfully. Syntax analysis heatmap is located at $OUTPUT_PATH/koala-stx-analysis.pdf"
echo '---------------------------------------------------------------------------------------------------------------------------------------'
popd
