# Koala

> If you find any *issues*, [open a GitHub issue](https://github.com/binpash/benchmarks/issues/new/choose).

Koala is a benchmark suite aimed at the characterization of performance-oriented research targeting the POSIX shell (i.e. executable on Dash, Bash, Zsh, and many other environments).
It combines a systematic collection of diverse shell programs collected from tasks found out in the wild, various real inputs to these programs.

The suite has been evaluated as part of the [ATC'25 Artifact Evaluation process](https://www.usenix.org/conference/atc25/call-for-artifacts) and
has received all three badges (*Available*, *Functional*, and *Reproduced*). That version (frozen in time), can be found at the [atc25-ae branch](https://github.com/kbensh/koala/tree/atc25-ae).

If any aspect of the suite is useful, please use the following citation:
```bibtex
@inproceedings{koala2025atc,
  title = {The Koala Benchmarks for the Shell: Characterization and Implications},
  author = {Evangelos Lamprou and Ethan Williams and Georgios Kaoukis and Zhuoxuan Zhang
        and Michael Greenberg and Konstantinos Kallas and Lukas Lazarek and Nikos Vasilakis},
  booktitle = {Proceedings of the 2025 USENIX Annual Technical Conference (USENIX ATC '25)},
  year = {2025},
  address = {Santa Clara, CA},
  publisher = {USENIX Association},
}
```

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

  - [Benchmarks](#benchmarks)
  - [Instructions](#instructions)
  - [Philosophy](#philosophy)
  - [Using the Suite](#using-the-suite)
    - [Environment & Setup Notes](#environment--setup-notes)
    - [Configuration](#configuration)
    - [Files produced per run](#files-produced-per-run)
    - [Usage examples](#usage-examples)
    - [Dynamic Characterization & Analysis](#dynamic-characterization--analysis)
      - [Running the Dynamic Analysis Separately](#running-the-dynamic-analysis-separately)
      - [Step-by-step:](#step-by-step)
      - [Anatomy of stats file](#anatomy-of-stats-file)
    - [Static Characterization & Analysis](#static-characterization--analysis)
  - [Contributing](#contributing)
- [Appendix I: Inputs](#appendix-i-inputs)
- [Appendix II: Dependencies](#appendix-ii-dependencies)
- [License](#license)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Benchmarks

| Benchmark   | Description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| `analytics`   | processes real-world network logs to extract and summarize key events.          |
| `bio`         | performs genomic and transcriptomic analysis using population and RNA-seq data. |
| `ci-cd`       | builds and tests open-source software projects.                             |
| `covid`       | analyzes public transit activity during the covid-19 pandemic.              |
| `file-mod`    | compresses, encrypts, and converts various file formats.                    |
| `inference`   | runs media-related inference tasks using large foundation models.           |
| `ml`          | implements a full machine learning pipeline using scikit-learn.             |
| `nlp`         | processes books using shell-based nlp pipelines from unix for poets.        |
| `oneliners`   | executes classic and modern one-liner shell pipelines.                      |
| `pkg`         | builds aur packages and analyzes npm packages for permissions.              |
| `repl`        | performs security auditing and git-based development workflow replay.       |
| `unixfun`     | solves unix text-processing problems from the 50-year anniversary challenge.|
| `weather`     | computes and visualizes historical weather statistics.                      |
| `web-search`  | implements crawling, indexing, and querying of wikipedia data.              |

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

## Using the Suite

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

# For development, mount the koala directory
$ sudo docker run -it --rm --cap-add=SYS_PTRACE --cap-add=NET_RAW --cap-add=NET_ADMIN --security-opt seccomp=unconfined --security-opt apparmor=unconfined -v "$(pwd):/koala" ghcr.io/binpash/benchmarks:latest bash
```

### Configuration
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

**Extra files produced when `-n <N>` > 1**
| File (aggregated)                     | Description                                             | Requires flag |
|---------------------------------------|---------------------------------------------------------|---------------|
| `<prefix>_stats_aggregated.txt`       | Mean / min / max of every numeric resource metric.      | `--resources` |
| `<benchmark>_times_aggregated.txt`    | Mean / min / max of wall-clock seconds.                 | `--time`      |

### Usage examples

1. Plain correctness run for the `unixfun` benchmark:
```bash
./main.sh unixfun
```
2. Run 10 times, record runtimes only:
```bash
./main.sh unixfun -n 10 --time
```
3. Heavy resource tracing inside Docker – 3 repetitions:
```bash
./main.sh unixfun -n 3 --resources
```
4. Lightweight local resource logging (no Docker):
```bash
./main.sh unixfun --resources --bare
```
5. Combine timing + resources, forward extra args to benchmark's infrastructure scripts:
```bash
./main.sh unixfun -n 5 --resources --time -- --small --fast
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
==================================================
```

Per-input-byte numbers are computed automatically: if `BENCHMARK_INPUT_FILE` points to a file or a directory, the harness figures out its byte size.
For more accurate analysis, please run inside a docker container.

Local (`--bare`) and Docker-based stats share the exact same format, so they aggregate seamlessly.

### Static Characterization & Analysis

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

## Contributing

Contributions are always welcome! Further details on how to contribute to the Koala benchmark suite project,
take a look at the [CONTRIBUTING](./CONTRIBUTING.md) file.

# Appendix I: Inputs

The table below contains all links to the inputs. Note: Some of these inputs are extremely large and hosted on low-bandwidth permanent-storage services such as Zenodo.

| Benchmark       | Min                                                                                                   | Small                                                                                                             | Full                                                                                                                                     |
| --------------- | ----------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `analytics`     | [Github Repo](https://github.com/binpash/benchmarks/tree/main/log-analysis/min_inputs)                | [Brown](https://atlas.cs.brown.edu/data/pcaps.zip) [Zenodo](https://zenodo.org/records/15361083)                  | [Brown](https://atlas.cs.brown.edu/data/pcaps_large.zip) [Zenodo](https://zenodo.org/records/15368510)                                   |
| `bio`           | [Github Repo](https://github.com/binpash/benchmarks/tree/main/bio/min_inputs)                         | [Brown](https://atlas.cs.brown.edu/data/bio/small) [Zenodo](https://zenodo.org/records/15361083)                  | [Brown](https://atlas.cs.brown.edu/data/bio/full) [Zenodo](https://zenodo.org/records/15367723)                                          |
| `ci-cd`         | No inputs                                                                                             | No inputs                                                                                                         | No inputs                                                                                                                                |
| `covid`         | [Github Repo](https://github.com/binpash/benchmarks/tree/main/covid-mts/min_inputs)                   | [Brown](https://atlas.cs.brown.edu/data/covid-mts/in_small.csv.gz) [Zenodo](https://zenodo.org/records/15361083)  | [Brown](https://atlas.cs.brown.edu/data/covid-mts/in_full.csv.gz) [Zenodo](https://zenodo.org/records/15368074)                          |
| `file-mod`      | [Github Repo](https://github.com/binpash/benchmarks/tree/main/file-enc/min_inputs)                    | [Brown](https://atlas.cs.brown.edu/data/pcaps.zip) [Zenodo](https://zenodo.org/records/15361083)                  | [Brown](https://atlas.cs.brown.edu/data/pcaps_large.zip) [Zenodo](https://zenodo.org/records/15368510)                                   |
| `inference`     | [Brown](https://atlas.cs.brown.edu/data/inference)                                                    | [Brown](https://atlas.cs.brown.edu/data/inference)                                                                                                              | [Brown](https://atlas.cs.brown.edu/data/inference)                                                                                                                                     |
| `weather`       | [Github Repo](https://github.com/binpash/benchmarks/tree/main/max-temp/min_inputs)                    | [Brown](https://atlas.cs.brown.edu/data/max-temp/noaa/) [Zenodo](https://zenodo.org/records/15361083)             | [Brown](https://atlas.cs.brown.edu/data/max-temp/noaa/) [Zenodo](https://zenodo.org/records/15368510)                                    |
| `ml`            | [Github Repo](https://github.com/binpash/benchmarks/tree/main/sklearn/inputs/covertype)               | [Brown](https://atlas.cs.brown.edu/data/sklearn/) [Zenodo](https://zenodo.org/records/15361083)                   | [Brown](https://atlas.cs.brown.edu/data/sklearn/) [Zenodo](https://zenodo.org/records/15368512)                                          |
| `nlp`           | [Brown](https://atlas.cs.brown.edu/data/gutenberg)                                                    | [Brown](https://atlas.cs.brown.edu/data/gutenberg/) [Zenodo](https://zenodo.org/records/15361083)                 | [Brown](https://atlas.cs.brown.edu/data/gutenberg/) [Zenodo](https://zenodo.org/records/15368510)                                        |
| `oneliners`     | [Brown](https://atlas.cs.brown.edu/data/dummy/)                                                       | [Brown](https://atlas.cs.brown.edu/data/dummy/) [Zenodo](https://zenodo.org/records/15361083)                     | [Brown](https://atlas.cs.brown.edu/data/dummy/) [Zenodo](https://zenodo.org/records/15368512)                                            |
| `pkg`           | [Github Repo](https://atlas.cs.brown.edu/data/packages_min)                                           | [Brown](https://atlas.cs.brown.edu/data/packages) [Zenodo](https://zenodo.org/records/15361083)                   | [Brown](https://atlas.cs.brown.edu/data/packages) [Zenodo](https://zenodo.org/records/15367723)                                          |
| `repl`          | No inputs                                                                                             | No inputs                                                                                                         | No inputs                                                                                                                                |
| `unixfun`       | [Github Repo](https://atlas.cs.brown.edu/data/unixfun)                                                 | [Brown](https://atlas.cs.brown.edu/data/unixfun/small) [Zenodo](https://zenodo.org/records/15361083)               | [Brown](https://atlas.cs.brown.edu/data/unixfun/large) [Zenodo](https://zenodo.org/records/15368512)                                      |
| `web-search`    | [Brown](https://atlas.cs.brown.edu/data/wikipedia_min.tar.gz)                                         | [Brown](https://atlas.cs.brown.edu/data/wikipedia_small.tar.gz) [Zenodo](https://zenodo.org/records/15361083)     | [Brown](https://atlas.cs.brown.edu/data/wikipedia.tar.gz) [Zenodo](https://zenodo.org/records/15368512)                                  |

# Appendix II: Dependencies 

Each benchmark includes dependencies across three categories: (1) software packages, (2) input datasets, and (3) miscellaneous dependencies (eg., endpoints, keys, etc.).
All inputs are permanently stored and available online.

| Benchmark        | Type       | Dependencies                                                                                                                                                                                                                                                                                   |
| --------------   | ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                                                                                 |
| `analytics`      | Inputs     | Log data (.pcap, .pcapng, and text files)                                                                                                                                                                                                                                                      |
| `analytics`      | Packages   | tcpdump, minimap2, bcftools, python3-pip, vim, ffmpeg, unrtf, imagemagick, zstd, liblzma-dev, libbz2-dev, zip, nodejs                                                                                                                                                                          |
| `bio`            | Inputs     | Genome inputs (.bam files)                                                                                                                                                                                                                                                                     |
| `bio`            | Packages   | samtools                                                                                                                                                                                                                                                                                       |
| `ci-cd`          | Inputs     | makeself git repository  Git repositories (memcached, redis, sqlite, vim, xz)                                                                                                                                                                                                                  |
| `ci-cd`          | Packages   | binutils, build-essential, coreutils, make, pbzip2, binutils, bzip2, zstd, gnupg, autotools-dev, automake, build-essential, clang, gcc, libevent-dev, libice-dev, libreadline-dev, libselinux-dev, libsm-dev, libtool, libtool-bin, libx11-dev, libxdmcp-dev, libxt-dev, make, pkg-config, tcl |
| `covid`          | Inputs     | Schedule data (.csv files)                                                                                                                                                                                                                                                                     |
| `covid`          | Packages   | python3-pip, libarchive-tools, libncurses5-dev, libncursesw5-dev, zstd, liblzma-dev, libbz2-dev, zip                                                                                                                                                                                           |
| `file-mod`       | Inputs     | Target files (.pcap and .pcapng files) Media files (.jpg and .wav files)                                                                                                                                                                                                                       |
| `file-mod`       | Packages   | apt-get, ffmpeg, unrtf, imagemagick, libarchive-tools, libncurses5-dev, libncursesw5-dev, zstd, liblzma-dev, libbz2-dev, zip, nodejs, tcpdump, ffmpeg, unrtf, imagemagick, libarchive-tools, libncurses5-dev, libncursesw5-dev, zstd, liblzma-dev, libbz2-dev, zip, xz-utils                   |
| `inference`      | Inputs     | Media files (.wav, .jpg) Model files                                                                                                                                                                                                                                                                                           |
| `ml`             | Inputs     | Text data (.csv files)                                                                                                                                                                                                                                                                         |
| `ml`             | Packages   | python3, python3-pip, joblib==1.4.2, numpy==1.26.4, scikit-learn==1.5.0, scipy==1.13.1, threadpoolctl==3.5.0                                                                                                                                                                                   |
| `nlp`            | Inputs     | Text data (.txt files)                                                                                                                                                                                                                                                                         |
| `oneliners`      | Inputs     | Text data (.txt files)                                                                                                                                                                                                                                                                         |
| `oneliners`      | Packages   | bsdmainutils, file, dos2unix                                                                                                                                                                                                                                                                   |
| `pkg`            | Inputs     | 150 AUR PKGBUILDs                                                                                                                                                                                                                                                                              |
| `pkg`            | Misc       | makedeb public key                                                                                                                                                                                                                                                                             |
| `pkg`            | Packages   | gpg, ffmpeg, unrtf, imagemagick, libarchive-tools, libncurses5-dev, libncursesw5-dev, zstd, liblzma-dev, libbz2-dev, makedeb, curl                                                                                                                                                             |
| `repl`           | Misc       | https://api.ipify.org endpoint                                                                                                                                                                                                                                                                 |
| `repl`           | Packages   | bash, fail2ban, gawk, grep, iproute2, iptables, net-tools, procps, ufw                                                                                                                                                                                                                         |
| `unixfun`        | Inputs     | Text data (.txt files)                                                                                                                                                                                                                                                                         |
| `weather`        | Inputs     | Temperature data (ISD files)                                                                                                                                                                                                                                                                   |
| `weather`        | Packages   | python3-pip, unrtf, libarchive-tools, libncurses5-dev, libncursesw5-dev, zstd, liblzma-dev, libbz2-dev, zip                                                                                                                                                                                    |
| `web-search`     | Inputs     | Wikipedia dump (.html files), Index (.txt file)                                                                                                                                                                                                                                                |
| `web-search`     | Packages   | p7zip-full, nodejs, npm, pandoc, natural                                                                                                                                                                                                                                                       |


# License

The Koala Benchmarks are licensed under the MIT License. See the LICENSE file for more information.
