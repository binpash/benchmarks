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

### Docker
```sh
# Build the container
$ docker build -t koala .

# Run the container
$ docker run --cap-add NET_ADMIN --cap-add NET_RAW -it koala

# For development, mount the benchmarks directory
$ docker run --cap-add NET_ADMIN --cap-add NET_RAW -it -v "$(pwd):/benchmarks" koala
```

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

### Dynamic Characterization & Analysis
When you run a benchmark with the `--resources` flag, existing process logs in  
`infrastructure/target/process-logs/` are automatically moved to  
`infrastructure/target/backup-process-logs/`.  
This ensures clean output when collecting new resource statistics. The new logs are then used to generate dynamic analysis visualizations.

#### Running the Dynamic Analysis Separately
You can also **run the dynamic analysis independently** of the main harness.  
This is useful for manually generating plots and debugging, while only having to execute each benchmark and avoiding the full setup and cleanup process.

#### Step-by-step:
1. **Install dependencies**:
   ```bash
   sudo apt-get install -y autoconf automake libtool build-essential cloc
   ```
   `pip install --break-system-packages -r "infrastructure/requirements.txt"`

2. **Run the analysis script manually:**:
    `./infrastructure/run_dynamic.py benchmark_name`
    This generates new process logs in:
    `infrastructure/target/process-logs/`

3. **(If using Docker):**  
   Copy the log files from the container to your host system in order to generate plots. The requirements will need to be installed in your host machine as well.  
   If you'd prefer not to copy files, you can instead run the visualizer with the `--text` flag to produce textual output directly:

   ```bash
   infrastructure/viz/dynamic.py --text
   ```

4. **Navigate to the infrastructure directory**:
    `cd infrastructure`

5. **Delete previous analysis output**:
    `rm -f target/dynamic_analysis.csv`

6. **Regenerate the analysis CSV**:
    `make target/dynamic_analysis.csv`

7. **Generate the visualizations**:
    `python infrastructure/viz/dynamic.py /path/to/output`

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
```
Per-input-byte numbers are computed automatically: if `BENCHMARK_INPUT_FILE` points to a file or a directory, the harness figures out its byte size. For more accurate analysis, please run inside a docker container.

Local (`--bare`) and Docker-based stats share the exact same format, so they aggregate seamlessly.

### Syntactic Characterization & Analysis

We use [`libdash`](https://github.com/binpash/libdash) to parse and analyze both the shell portion of each benchmark, and the portions of components called into by the shell and which often implement the kernel of a computation: for the shell portion, we count the total occurrences of every AST node; for the command portion, we analyze only AST nodes counting commands, built-ins, and functions—noting that the results are conservative, as they do not count dynamic commands.

The analysis produces CSV summaries and heatmaps across the benchmark suite, highlighting the use of each shell construct.

1. **Install dependencies**:
   ```bash
   sudo apt-get install -y autoconf automake libtool build-essential cloc
   pip install --break-system-packages -r infrastructure/requirements.txt

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