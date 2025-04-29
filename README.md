# Koala

> _It's a suite that has benchmarks in it._

## Benchmarks

| Benchmark    | Description                                             |
| ---------    | -----------                                             |
| aurpkg       |                                                         |
| bio          | Bioinformatics.                                         |
| covid-mts    | COVID-19 multivariate time series.                      |
| dpt          | Hieroglyph classification benchmark.                    |
| file-enc     | File encoding.                                          |
| git-workflow | A Git workflow.                                         |
| llm          | Various benchmarks using llms.                          |
| log-analysis | Log analysis.                                           |
| makeself     | Make self-extractable archives on Unix                  |
| max-temp     | Maximum temperature.                                    |
| media-conv   | Media conversion.                                       |
| nlp          | Natural language processing.                            |
| oneliners    | One-liners.                                             |
| opt-parallel | Fast map-reduce equivalent.                             |
| port-scan    | Port-scan processing.                                   |
| ray-tracing  | Ray-tracing log processing.                             |
| riker        |                                                         |
| sklearn      | Machine learning.                                       |
| teraseq      |                                                         |
| tuft-weather | Tuft weather.                                           |
| uniq-ips     | Unique IPs.                                             |
| unix50       | Unix 50.                                                |
| web-index    | Web index.                                              |

## Instructions
First, set the shell runtime to benchmark with the `$BENCHMARK_SHELL` shell variable.
By default, this is set to `bash`.
You can also pass in flags/options to use with the shell.
For example, to benchmark PaSh with `--width 4`, run `export BENCHMARK_SHELL="$PASH_TOP/pa.sh --width 4"`.

`main.sh` is the one-stop harness for downloading dependencies and inputs, running, profiling and verifying a **single benchmark** in this suite.

```bash
./main.sh <BENCHMARK_NAME> [OPTIONS] [-- args passed to run.sh …]
```

`run_all.sh` runs the `main.sh` script for all benchmarks inside a Docker container or locally if the `--bare` flag is used.

### Core options

| Flag / Option                         | Effect                                                                                                              | Typical use-case                                   |
|---------------------------------------|---------------------------------------------------------------------------------------------------------------------|----------------------------------------------------|
| **`-n <N>` / `--runs <N>`**           | Execute the benchmark **N** times (default = 1).                                                                    | Measure variance, warm-up caches, find flaky runs. |
| **`--resources`**                     | Collect CPU/RAM/I/O stats.<br>Writes `*_stats_run<i>.txt` for every run **and** a summary `*_stats_aggregated.txt`. | Profiling, optimisation, capacity planning.        |
| **`--bare`**                          | Use the lightweight local logger instead of the Docker-based tracer.                                                | When Docker isn’t available or is too heavy.       |
| **`-t` / `--time`**                   | Measure wall-clock runtime with `/usr/bin/time`.<br>Produces `{benchmark}_times_aggregated.txt` when `-n > 1`.      | Quick speed checks, perf regression testing.       |
| **`--small`**                         | Run the benchmark with a reduced (small) input set.                                                                 | Fast experiments, small-scale characterisation.    |
| **`--min`**                           | Run the benchmark with the absolute minimum inputs.                                                                 | Suite-level sanity checks.                         |

Flags, apart from those referring to input sizes, can be combined freely (e.g. `--resources --bare -n 5`).

### Files produced per run

| File (per-run)                                | Contents / Purpose                                                                                 | Generated when …                                              |
|-----------------------------------------------|----------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| `<benchmark>.out` / `<benchmark>.err`         | Stdout / stderr from `run.sh`.                                                                     | **Always**                                                    |
| `benchmark.hash`                              | Pass/fail hashes written by `verify.sh` – indicates whether the run was correct.                   | **Always**                                                    |
| `logs/…` (folder)                             | Raw monitors: `*.pidstat`, `*.io`, `*.cpu`, `*.mem`, `*.time`, `*.val`.                            | Only with (**`--resources`** or **`--time`**) &&  **`--bare`**|
| `<prefix>_stats_run<i>.txt`                   | Human-readable CPU/RAM/I/O summary for run *i*.                                                    | Only with **`--resources`**                                   |
| `<benchmark>_time_run<i>.val`                 | Single wall-clock number (seconds) for run *i*.                                                    | Only with **`--time`**                                        |

### Extra files produced when **`-n <N>`** > 1

| File (aggregated)                     | Description                                             | Requires flag |
|---------------------------------------|---------------------------------------------------------|---------------|
| `<prefix>_stats_aggregated.txt`       | Mean / min / max of every numeric resource metric.      | `--resources` |
| `<benchmark>_times_aggregated.txt`    | Mean / min / max of wall-clock seconds.                 | `--time`      |

### Usage examples

#### 1. Plain correctness run (old behaviour)
```bash
./main.sh unix50
```
#### 2. Run 10×, record runtimes only
```bash
./main.sh unix50 -n 10 --time
```
#### 3. Heavy resource tracing inside Docker – 3 repetitions
```bash
./main.sh unix50 -n 3 --resources
```
#### 4. Lightweight local resource logging (no Docker)
```bash
./main.sh unix50 --resources --bare
```
#### 5. Combine timing + resources, forward extra args to benchmark's infrastructure scripts
```bash
./main.sh unix50 -n 5 --resources --time -- --small --fast
```
### Anatomy of stats file
```
Benchmark Statistics
==================================================
Benchmark: bio
--------------------------------------------------
Total CPU time: 6.68 sec
Total Wall time: 6.91 sec
Total IO bytes: 1107624733.00
Max Memory Usage: 17235968.00 bytes
CPU time per input byte: 0.000000 sec/byte
Memory per input byte: 0.071815 bytes/byte
IO per input byte: 4.615032 bytes/byte
Time in Shell: 0.00 sec
Time in Commands: 6.68 sec
```
Per-input-byte numbers are computed automatically: if `BENCHMARK_INPUT_FILE` points to a file or a directory, the harness figures out its byte size. For more accurate analysis, please run inside a docker container.

Local (`--bare`) and Docker-based stats share the exact same format, so they aggregate seamlessly.

### Docker
```sh
# Build the container
$ docker build -t koala .

# Run the container
$ docker run --cap-add NET_ADMIN --cap-add NET_RAW -it koala

# For development, mount the benchmarks directory
$ docker run --cap-add NET_ADMIN --cap-add NET_RAW -it -v "$(pwd):/benchmarks" koala
```
