# Transit Analytics Benchmarks

This directory contains the Transit Analytics benchmarks.

## TLDR
```bash
./setup.sh
python3 run --window 10 --target both
```

Results are in `<hs_top>/report/output/bus-analytics`.

## Directory Structure

The benchmark directory is structured with individual subdirectories for each benchmark, identified by numbers. Each subdirectory contains a shell script named after the benchmark number (e.g., `1.sh` for benchmark 1) and a `run` script to execute the benchmark.

```
.
├── 1
│   ├── 1.sh
│   └── run
├── 2
│   ├── 2.sh
│   └── run
...
├── full.sh
├── run
└── setup.sh
```

## Setup

Before running the benchmarks, you need to set up the input data. This setup is performed by the `setup.sh` script located in the root of the benchmark directory. To run the setup, simply execute:

```sh
./setup.sh
```

This script prepares the necessary input data for all benchmarks.

The input files are downloaded and stored in the `<hs_top>/report/resources/bus-analytics` directory. Make sure you have sufficient space and permissions in this directory before running the setup.


## Running Benchmarks

To execute the entire suite of benchmarks, run the `run` script located in the root of the benchmark directory with Python 3:

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

## Individual Benchmarks

You can also run specific benchmarks individually by navigating to the subdirectory for the benchmark and running the `run` script with the same options. For example, to run benchmark number 1:

```bash
cd 1
python3 run --window 10 --target hs-only --log enable
```

Replace `1` with the appropriate benchmark number as needed.


The outputs from running the benchmarks are saved in the `<hs_top>/report/output/bus-analytics` directory. Each sub-benchmark will have its own subdirectory. Each subdirectory will contain the following files:

- `error`: empty if the benchmark execution completes successfully without differences in output between `sh` and `hs`; non-empty if differences are found.
- `hs_time`: execution time of the `hs`.
- `sh_time`: execution time of the `sh`.
- `hs_log`: logs of `hs`.


### Benchmark details
 Mass-Transport System Analytics

This set of scripts script is part of [a recent study on OASA](https://insidestory.gr/article/noymera-leoforeia-athinas) from Diomidis Spinellis and Eleftheria Tsaliki.  OASA is the the mass-transport system supporting the city of Athens. 

1. `1.sh`: Vehicles on the road per day
2. `2.sh`: Days a vehicle is on the road
3. `3.sh`: Hours each vehicle is on the road
4. `4.sh`: Hours monitored each day
5. `5.sh`: Hours each bus is active each day


