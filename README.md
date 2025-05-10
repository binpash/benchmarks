# Koala

> Contact: `koala@brown.edu` or [open a GitHub issue](https://github.com/binpash/benchmarks/issues/new/choose).

Koala is a benchmark suite aimed at the characterization of performance-oriented research targeting the POSIX shell (i.e. executable on Dash, Bash, Zsh, and many other environments).
It combines a systematic collection of diverse shell programs collected from tasks found out in the wild, various real inputs to these programs.

> **Artifact evaluation:**
>
> To evaluate key results in our ATC'25 paper titled _"The Koala Benchmarks for the Shell: Characterization and Implications"_, jump straight to [`INSTRUCTIONS/`](https://github.com/binpash/benchmarks/blob/main/INSTRUCTIONS.md).

## Benchmarks

| Benchmark    | Description                                             |
| ---------    | -----------                                             |
| `aurpkg`       | AUR package builds.                                     |
| `bio`          | Bioinformatics.                                         |
| `covid-mts`    | COVID-19 multivariate time series.                      |
| `file-enc`     | File encoding.                                          |
| `log-analysis` | Log analysis.                                           |
| `makeself`     | Make self-extractable archives on Unix.                 |
| `max-temp`     | Maximum temperature.                                    |
| `media-conv`   | Media conversion.                                       |
| `nlp`          | Natural language processing.                            |
| `oneliners`    | One-liners.                                             |
| `riker`        | Incremental builds.                                     |
| `sklearn`      | Machine learning.                                       |
| `vps-audit`    | Audit a Linux machine.                                  |
| `unix50`       | Unix 50.                                                |
| `web-index`    | Web index.                                              |

## Instructions

The top-level `main.sh` script is a quick script for downloading dependencies and inputs, running, profiling, and verifying a _single Koala benchmark_.

```bash
./main.sh <BENCHMARK_NAME> [OPTIONS] [<args passed to execute.sh>]
```

## Philosophy

To support the diverse landscape of shell programs and shell-related research, the Koala benchmark suite is designed with flexibility and simplicity in mind.

Its infrastructure is deliberately minimal and easy to modify, making it adaptable to a wide variety of systems and use cases.
As it cannot anticipate all potential applications, Koala encourages users to modify any part of the infrastructure to better suit their needs.

For example, given the following (example) benchmark:

```sh
# example-benchmark/scripts/x.sh
cat file.txt | grep "foo" | wc -l
```

Someone experimenting with GNU parallel might want to modify the script as follows:

```sh
# example-benchmark/scripts/x.sh
cat file.txt | parallel --pipe grep "foo" | wc -l
```

Similalary, a research team developing a [distributed shell](https://www.usenix.org/conference/nsdi23/presentation/mustafa) can modify the script to:

```sh
# example-benchmark/scripts/x.sh
hdfs dfs -cat file.txt | dsh --pipe grep "foo" | wc -l
```

For systems that act as a drop-in replacement for the shell can  use Koala's benchmarks by overriding the `$KOALA_SHELL` variable to point to their system.
For example, to apply [the PaSh system](https://www.usenix.org/conference/osdi22/presentation/kallas) to the Koala benchmarks, one can do:

```sh
$ export KOALA_SHELL="./pa.sh --width 4"
$ ./main.sh example-benchmark
```

### Benchmark Structure

Each benchmark directory contains:

```
benchmarks/<name>/
├── scripts/        # Benchmark scripts (*.sh)
├── install.sh      # Installs dependencies
├── fetch.sh        # Fetches input data
├── execute.sh      # Runs benchmark
├── validate.sh     # Validates output via hashes
└── clean.sh        # Cleans temporary files (input and output files)
```

To manually execute a single benchmark: 

```sh
cd benchmarks/<name>

./install.sh

# This will place input data in the `benchmarks/<name>/inputs` folder
./fetch.sh

# This will run one by one all scripts inside `benchmarks/<name>/scripts/` folder
# Any output files produced will be placed in the `benchmarks/<name>/outputs` folder
./execute.sh

# This will check the output files against the expected hashes, and print the results
# For benchmarks that do not produce output files, specialized validation logic is used
./validate.sh

# This will remove all temporary files created by the benchmark
# By default, it will remove both inputs and outputs, returning the benchmark folder in its original state
./clean.sh
```

The harness which automates this process and collects/displays metrics is `main.sh`:

```
benchmarks/main.sh              # Drives all benchmarks
```

Configuration options (via env or CLI):
To control the shell interpreter used to run the benchmarks, you can set the
`KOALA_SHELL` environment variable.

### Environment & Setup Notes

**Note:** The setup scripts in this suite are designed for _Debian-based systems_.
Koala comes with a Docker image, highly recommended when working on non-Debian systems.

To build and run the Docker image:
```sh
# Build the container
$ docker build -t koala .

# Run the container
$ docker run -it koala

# For development, mount the benchmarks directory
$ docker run -it -v "$(pwd):/benchmarks" koala
```

### Usage
```
Usage: ./main.sh BENCHMARK_NAME [--time|--resources|--bare|args...]
  --min            Run the benchmark with minimal inputs (default)
  --small          Run the benchmark with small inputs
  --full          Run the benchmark with full inputs
  --time, -t       Measure wall-clock time
  --resources      Measure resource usage
  --bare           Run locally without Docker
  --runs, -n N     Number of runs (default: 1)
  --clean, -c      Run the full cleanup script (both inputs and outputs)
  --keep, -k       Keep outputs
  --prune          Run the benchmark on a fresh container (will need to re-download everything on each run)
  --help, -h       Show this help message
```

Flags, apart from those referring to input sizes, can be combined freely (e.g. `--resources --bare -n 5`).

### Files produced per run
| File (per-run)                                | Contents / Purpose                                                                                 | Generated when …                                              |
|-----------------------------------------------|----------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| `<benchmark>.out` / `<benchmark>.err`         | Stdout / stderr from `execute.sh`.                                                                 | **Always**                                                    |
| `benchmark.hash`                              | Pass/fail hashes written by `validate.sh` – indicates whether the run was correct.                 | **Always**                                                    |
| `logs/…` (folder)                             | Raw monitors: `*.pidstat`, `*.io`, `*.cpu`, `*.mem`, `*.time`, `*.val`.                            | Only with (**`--resources`** or **`--time`**) &&  **`--bare`**|
| `<prefix>_stats_run<i>.txt`                   | Human-readable CPU/RAM/I/O summary for run *i*.                                                    | Only with **`--resources`**                                   |
| `<benchmark>_time_run<i>.val`                 | Single wall-clock number (seconds) for run *i*.                                                    | Only with **`--time`**                                        |

### Extra files produced when **`-n <N>`** > 1
| File (aggregated)                     | Description                                             | Requires flag |
|---------------------------------------|---------------------------------------------------------|---------------|
| `<prefix>_stats_aggregated.txt`       | Mean / min / max of every numeric resource metric.      | `--resources` |
| `<benchmark>_times_aggregated.txt`    | Mean / min / max of wall-clock seconds.                 | `--time`      |

### Usage examples

1. Plain correctness run for the `unix50` benchmark:
```bash
./main.sh unix50
```
2. Run 10 times, record runtimes only:
```bash
./main.sh unix50 -n 10 --time
```
3. Heavy resource tracing inside Docker – 3 repetitions:
```bash
./main.sh unix50 -n 3 --resources
```
4. Lightweight local resource logging (no Docker):
```bash
./main.sh unix50 --resources --bare
```
5. Combine timing + resources, forward extra args to benchmark's infrastructure scripts:
```bash
./main.sh unix50 -n 5 --resources --time -- --small --fast
```

### Dynamic Characterization & Analysis
When running a benchmark with the `--resources` flag, existing process logs in `infrastructure/target/process-logs/` are automatically moved to `infrastructure/target/backup-process-logs/`.
This ensures clean output when collecting new resource statistics. The new logs are then used to generate dynamic analysis visualizations.

#### Running the Dynamic Analysis Separately
You can also run the dynamic analysis independently of the main harness. This is useful for manually generating plots and debugging, while only having to execute each benchmark and avoiding the full setup and cleanup process.

#### Step-by-step:
1. **Install dependencies**:
   ```bash
   sudo apt-get install -y autoconf automake libtool build-essential cloc
   pip install --break-system-packages -r "infrastructure/requirements.txt"
    ```

2. **Run the analysis script manually**:
    ```bash
    ./infrastructure/run_dynamic.py benchmark_name
    ```
    This generates new process logs in:
    ```bash
    infrastructure/target/process-logs/
    ```

3. **(If using Docker)**:  
   Copy the log files from the container to your host system in order to generate plots. The requirements will need to be installed in your host machine as well.  
   If you'd prefer not to copy files, you can instead run the visualizer with the `--text` flag to produce textual output directly:

   ```bash
   infrastructure/viz/dynamic.py --text
   ```

4. **Navigate to the infrastructure directory**:
    ```bash
    cd infrastructure

5. **Delete previous analysis output**:
    ```bash
    rm -f target/dynamic_analysis.csv
    ```

6. **Regenerate the analysis CSV**:
    ```bash
    make target/dynamic_analysis.csv
    ```

7. **Generate the visualizations**:
    ```bash
    python infrastructure/viz/dynamic.py /path/to/output
    ```

This produces benchmark-specific performance plots, showing shell vs command time,
CPU usage, I/O throughput, and memory footprint, for all benchmarks that have logs present in `infrastructure/target/process-logs/`

#### Anatomy of stats file
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
--------------------------------------------------
```

Per-input-byte numbers are computed automatically: if `BENCHMARK_INPUT_FILE` points to a file or a directory, the harness figures out its byte size.
For more accurate analysis, please run inside a docker container.

Local (`--bare`) and Docker-based stats share the exact same format, so they aggregate seamlessly.

### Syntactic Characterization & Analysis

We use [`libdash`](https://github.com/binpash/libdash) to parse and analyze both the shell portion of each benchmark, and the portions of components called into by the shell and which often implement the kernel of a computation: for the shell portion, we count the total occurrences of every AST node; for the command portion, we analyze only AST nodes counting commands, built-ins, and functions—noting that the results are conservative, as they do not count dynamic commands.

The analysis produces CSV summaries and heatmaps across the benchmark suite, highlighting the use of each shell construct.

1. **Install dependencies**:
   ```bash
   sudo apt-get install -y autoconf automake libtool build-essential cloc
   pip install --break-system-packages -r infrastructure/requirements.txt
   ```
   
2. **Register the benchmark script**:  
   Add the new benchmark’s script pattern to:
   `infrastructure/data/script-globs.json`
   > **Note:** Syntactic analysis only works for **POSIX-compliant** scripts.

3. **Remove previous analysis artifacts**:
    ```bash
    rm -f infrastructure/target/cyclomatic.csv
    rm -f infrastructure/target/lines_of_code.csv
    rm -f infrastructure/target/nodes_in_scripts.csv
    rm -f infrastructure/target/scripts_to_benchmark.csv
    ```

4. **Navigate to the infrastructure directory**:
    `cd infrastructure`

5. **Regenerate the syntactic analysis artifacts**:
    `make`

6. **Generate visualizations**:
    ```
    python infrastructure/viz/syntax.py output_dir
    python infrastructure/viz/commands.py output_dir
    ```

These will produce plots summarizing shell syntax usage and external command invocation patterns for all registered benchmarks in the specified `output_dir`.

## License

The Koala Benchmarks are licensed under the MIT License. See the LICENSE file for more information.
