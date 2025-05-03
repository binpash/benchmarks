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

* Permanently hosted at: [https://github.com/blind/koala](https://github.com/blind/koala)
* MIT-licensed and open-source
* Includes containerized and bare-metal support

To download and execute the minimal version on your machine:

```bash
curl -s koala-blind | sh && ./main.sh --min --bare
```

# Artifact Functional

## Documentation

KOALA includes documentation and reproducibility support:

* Each benchmark follows a modular layout: `install.sh`, `fetch.sh`, `execute.sh`, `validate.sh`, `clean.sh`
* Top-level driver: `main.sh` orchestrates execution
* Container support via `Dockerfile`
* Benchmarks include realistic input sizes hosted on institutional servers and archival storage

## Completeness

At a high level, the paper claims the following contributions (p. 2):

1. A set of real-world shell programs across domains (CI, ML, SysAdmin, etc.)
2. Accompanying inputs under three tiers: `--min`, `--small`, and `--full`
3. Automation of dependency setup, input validation, and correctness checking
4. Infrastructure for static and dynamic analysis and characterization of the benchmark suite
5. Application of the benchmark suite to a set of prior shell optimization tools

## Exercisability

**Quickstart (minimal test):**

```bash
curl -s koala-blind | sh
cd benchmarks
./main.sh --min --bare
```

Expected time: ~2.5 minutes. Runs all benchmarks with synthetic inputs and checks output hashes.

**Benchmark-specific example (nlp):**

```bash
cd benchmarks
./main.sh --small nlp
```

**Configurable system integration:**

```bash
KOALA_SHELL=./your_system.sh ./main.sh --small max-temp
```

**Containerized run:**

```bash
docker build -t koala .
docker run -v $(pwd):/koala -it koala /koala/main.sh --small
```

# Results Reproducible

The artifact supports evaluation of:

* Correctness: `validate.sh` uses output hashes to ensure functional correctness
* Performance: timings and metrics collected via `/proc` and `psutil`

**Reproducing full results (optional):**

```bash
./main.sh         # full-size inputs
./main.sh --small # smaller scale inputs
```

Execution may take several hours depending on hardware. Reports are auto-generated.

# Infrastructure Overview

Each benchmark directory contains:

```
benchmarks/<name>/
├── scripts/        # Benchmark scripts (e.g., *.sh)
├── install.sh      # Installs dependencies
├── fetch.sh        # Fetches input data
├── execute.sh      # Runs benchmark
├── validate.sh     # Validates output via hashes
└── clean.sh        # Cleans temporary files
```

Top-level:

```
main.sh              # Drives all benchmarks
Dockerfile           # Container support
```

Configuration options (via env or CLI):

* `KOALA_SHELL`: shell interpreter (default: `sh`)
* `KOALA_INFO`: metrics to collect (e.g., time, memory)

# Contact

For questions or bug reports, please contact the authors or open an issue on GitHub.
