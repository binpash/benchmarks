# Getting Started Instructions

The structure of this document mirrors the [USENIX ATC '25 artifact evaluation process](https://www.usenix.org/conference/atc25/call-for-artifacts). At a glance:

* [ ] [Artifact available](#artifact-available): Hosted on GitHub with container support
* [ ] [Artifact functional](#artifact-functional): Modular design, complete infrastructure, and usability across environments
* [ ] [Results reproducible](#results-reproducible): Correctness checks, benchmark suite execution, and reporting

To "kick the tires" of this artifact:

* skim this README for structure and setup paths (~2 minutes)
* jump to [Exercisability](#exercisability) for a quick demo (~3 minutes)

# Artifact Available

The KOALA benchmark suite is:

* Permanently hosted at: [https://github.com/binpash/koala](https://github.com/binpash/koala)
* MIT-licensed and open-source
* Includes containerized and bare-metal support

To download and execute the minimal version on your machine:

```sh
curl -s koala-blind | sh && ./main.sh --min --bare
```

# Artifact Functional

## Documentation

KOALA includes documentation and reproducibility support:

* Each benchmark follows a modular layout: `install.sh`, `fetch.sh`, `execute.sh`, `validate.sh`, `clean.sh`
* Top-level driver: `main.sh` orchestrates execution
* Container support via `Dockerfile`
* Benchmarks include realistic input sizes hosted on institutional servers and archival storage

### Benchmark Structure

Each benchmark directory contains:

```
benchmarks/<name>/
├── scripts/        # Benchmark scripts (e.g., *.sh)
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

* `KOALA_SHELL`: shell interpreter (default: `sh`)
* `KOALA_INFO`: metrics to collect (e.g., time, memory)

### Benchmark Suite

The suite contains 15 benchmarks. 

## Completeness

At a high level, the paper claims the following contributions (p. 2):

1. A set of real-world shell programs that span a variety of domains
2. Accompanying inputs under three tiers: `--min`, `--small`, and `--full`
3. Automation of dependency setup, input validation, and correctness checking
4. Infrastructure for static and dynamic analysis and characterization of the benchmark suite
5. Application of the benchmark suite to a set of prior shell optimization tools

## Exercisability

**Quickstart:**

```sh
git clone https://github.com/binpash/benchmarks
cd benchmarks
chmod +x kick-tires.sh
./kick-tires.sh
```

**Quickstart w/ Docker:**

```sh
docker pull ghcr.io/binpash/benchmarks:latest
docker run -it --rm -v "$(pwd)":/mnt ghcr.io/binpash/benchmarks:latest bash
```

**Quickstart w/ Docker (build):**

```sh
git clone https://github.com/binpash/benchmarks
cd benchmarks
docker build -t benchmarks .
docker run -it --rm -v "$(pwd)":/benchmarks benchmarks bash
```

Expected time: ~2.5 minutes. Runs all benchmarks with synthetic inputs and checks output hashes.

**Benchmark-specific example (nlp):**

```sh
cd benchmarks
./main.sh --small nlp
```

**Configurable system integration:**

To run each of the benchmarks on a system that acts as a shell interpreter, you
can set the `KOALA_SHELL` environment variable.

```sh
KOALA_SHELL=./your_system.sh ./main.sh --small max-temp
```

**Containerized run:**

```sh
docker build -t koala .
docker run -v $(pwd):/koala -it koala /koala/main.sh --small
```

# Results Reproducible

The artifact supports evaluation of:

* Correctness: `validate.sh` uses output hashes to ensure functional correctness
* Performance: timings and metrics collected via `/proc` and `psutil`

**Reproducing full results (optional):**

```sh
./main.sh         # full-size inputs
./main.sh --small # smaller scale inputs
./main.sh --min   # minimal inputs
```

Execution may take several hours depending on hardware. Reports are auto-generated.

# Contact

For questions or bug reports, please contact the authors or open an issue on GitHub.
