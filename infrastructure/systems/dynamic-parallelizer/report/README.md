# hs Benchmark Directory README

This is the benchmark directory of `hs`. This directory contains essential tools and scripts for running benchmarks, analyzing logs, generating reports, and visualizing results.

## Overview

The benchmarking tool offers a comprehensive interface to evaluate the performance of `hs` against `bash`, supporting a range of features:

- Execution of benchmarks using both `bash` and `hs`.
- Comparison of outputs and performance metrics between `Bash` and `hs`.
- Generation of detailed reports including Gantt charts for in-depth analysis.
- Creation of CSV files for data analysis and bar charts for visual comparison.
- Optional verbose output for detailed execution logs.

## Environment Variables

Several environment variables are essential for the benchmarking tool's operation:

- `ORCH_TOP`: The top directory of the orchestrator system.
- `WORKING_DIR`: Directory for benchmarks and reports.
- `TEST_SCRIPT_DIR`: Directory containing benchmark scripts.
- `RESOURCE_DIR`: Directory for storing resources required by benchmarks.
- `PASH_TOP`: The directory of `pash`.
- `PASH_SPEC_TOP`: The top directory of `hs`.

## Command-Line Interface

The benchmark runner's command-line interface includes options for controlling the output and behavior:

- `--no-plots`: Disables the generation of plot visualizations.
- `--no-logs`: Prevents saving log files.
- `--csv-output`: Enables saving results in CSV format.
- `--verbose`: Enables verbose output, providing detailed logs of the benchmarking process.
- `--full-gantt`: Generate a full Gantt chart for each benchmark.

## Benchmark Configuration

Benchmarks are configured via `benchmark_config.json`, with each entry specifying:

- `name`: Benchmark name.
- `env`: Environment variables required by the benchmark.
- `pre_execution_script`: Commands for initial setup, like fetching data.
- `command`: The benchmark command or script.
- `orch_args`: Arguments for the `hs` system.

Example configuration:

```json
[
    {
        "name": "Dgsh 1.sh - 120M",
        "env": ["INPUT_FILE={RESOURCE_DIR}/in120M.xml"],
        "pre_execution_script": ["wget -nc -O in120M.xml http://aiweb.cs.washington.edu/research/projects/xmltk/xmldata/data/dblp/dblp.xml"],
        "command": "{TEST_SCRIPT_DIR}/dgsh/1.sh",
        "orch_args": "-d 2 --sandbox-killing"
    }
]
```

## Running Benchmarks

To execute benchmarks:

1. Navigate to the benchmark runner directory (e.g., `cd ./report`).
2. Run the benchmark runner with the desired arguments (e.g., `python3 main.py --csv-output`).

Results, including logs, plots, and CSV files, are saved in the `report_output` directory.

## Results Interpretation

Results encompass:

- Execution times for `bash` and `hs`.
- Comparative analysis of execution times.
- Validity checks of outputs.
- Execution logs and error messages in verbose mode.
- Gantt charts for timeline analysis and bar charts for speculative execution analysis.

## Contributions

Contributions to enhance or expand the benchmark suite are welcome. When adding new benchmarks, ensure to update `benchmark_config.json` accordingly.