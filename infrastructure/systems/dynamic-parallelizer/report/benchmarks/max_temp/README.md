Based on your request to create a README for the `max_temp` benchmark similar to the DGSH benchmarks README, here's a suggested README markdown for the `max_temp` benchmark directory:

```markdown
# Max Temp Benchmark

This directory contains the Max Temp benchmark for DGSH.

## TLDR
```bash
./setup.sh
python3 run --window 10 --target both
```

Results are in `<hs_top>/report/output/max_temp`.

## Directory Structure

The `max_temp` benchmark directory contains the main script, a setup script, and a run script:

```
.
├── max_temp.sh
├── run
└── setup.sh
```

## Setup

The setup for the Max Temp benchmark prepares the necessary input data. Execute the setup script located in the root of the benchmark directory to prepare the input data:

```sh
./setup.sh
```

The input files are downloaded and stored in the `<hs_top>/report/resources/max_temp` directory. Ensure sufficient space and permissions are available in this directory before running the setup.

## Running the Benchmark

To run the Max Temp benchmark, use the `run` script located in the root of the benchmark directory with Python 3:

```bash
python3 run [-h] [--window WINDOW] [--target {hs-only,sh-only,both}] [--log {enable,disable}]
```

This script supports several options to customize the benchmark execution:

- `--window WINDOW`: Specifies the window size to use with `hs`.
- `--target {hs-only,sh-only,both}`: Determines whether to run benchmarks using `sh` only, `hs` only, or both.
- `--log {enable,disable}`: Enables or disables logging for `hs`.

For example, to run all benchmarks with a window size of 10, targeting both `sh` and `hs`, and with logging enabled:

```bash
python3 run --window 10 --target both --log enable
```

## Execution Outputs

The outputs from running the benchmarks are saved in the `<hs_top>/report/output/dgsh` directory. Each sub-benchmark will have its own subdirectory. Each subdirectory will contain the following files:

- `error`: empty if the benchmark execution completes successfully without differences in output between `sh` and `hs`; non-empty if differences are found.
- `hs_time`: execution time of the `hs`.
- `sh_time`: execution time of the `sh`.
- `hs_log`: logs of `hs`.
