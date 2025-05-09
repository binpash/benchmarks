# Overview

The structure of this document mirrors the [USENIX ATC '25 artifact evaluation process](https://www.usenix.org/conference/atc25/call-for-artifacts).
This artifact targets the following badges:

* [ ] [Artifact available](#artifact-available): Hosted on GitHub with container support
* [ ] [Artifact functional](#artifact-functional): Modular design, complete infrastructure, and usability across environments
* [ ] [Results reproducible](#results-reproducible): Correctness checks, benchmark suite execution, and reporting

To "kick the tires" of this artifact:

* skim this README for structure and setup paths (~2 minutes)
* jump to [Exercisability](#exercisability) for a quick demo (~20 minutes)

# Artifact Available

The KOALA benchmark suite is:

* Permanently hosted at: [https://github.com/binpash/koala](https://github.com/binpash/koala)

# Artifact Functional

## Documentation

[Repo README](https://github.com/binpash/benchmarks)

Per-benchmark documentation is available within each benchmark directory, it
contains information about the computation, input data, and expected output.
- [aurpkg](https://github.com/binpash/benchmarks/tree/main/aurpkg)
- [bio](https://github.com/binpash/benchmarks/tree/main/bio)
- [covid-mts](https://github.com/binpash/benchmarks/tree/main/covid-mts)
- [file-enc](https://github.com/binpash/benchmarks/tree/main/file-enc)
- [log-analysis](https://github.com/binpash/benchmarks/tree/main/log-analysis)
- [makeself](https://github.com/binpash/benchmarks/tree/main/makeself)
- [max-temp](https://github.com/binpash/benchmarks/tree/main/max-temp)
- [media-conv](https://github.com/binpash/benchmarks/tree/main/media-conv)
- [nlp](https://github.com/binpash/benchmarks/tree/main/nlp)
- [oneliners](https://github.com/binpash/benchmarks/tree/main/oneliners)
- [riker](https://github.com/binpash/benchmarks/tree/main/riker)
- [sklearn](https://github.com/binpash/benchmarks/tree/main/sklearn)
- [unix50](https://github.com/binpash/benchmarks/tree/main/unix50)
- [vps-audit](https://github.com/binpash/benchmarks/tree/main/vps-audit)
- [web-index](https://github.com/binpash/benchmarks/tree/main/web-index)

* Each benchmark follows a modular layout: `install.sh`, `fetch.sh`, `execute.sh`, `validate.sh`, `clean.sh`
* Top-level driver: `main.sh` orchestrates execution
* Container support via `Dockerfile`
* Benchmarks include realistic input sizes hosted on institutional servers and archival storage

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

To run the benchmarks directly on your machine, without containerization, use the following lines.
Note that the scripts *assume* some basic dependencies that very minimal system setups may not already have (e.g. `git`).

```sh
git clone https://github.com/binpash/benchmarks
cd benchmarks
sudo docker pull ghcr.io/binpash/benchmarks:latest
sudo docker tag ghcr.io/binpash/benchmarks:latest koala
sudo docker run -it --rm -v "$(pwd)":/benchmarks koala /bin/bash -c "cd /benchmarks && ./kick-tires.sh"
```

Expected time: ~2 -- 20 minutes, depending on hardware.
Runs all benchmarks with minimal inputs.
Most of the time will be spent on downloading inputs.

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
docker run -v $(pwd):/benchmarks -it koala /benchmarks/main.sh --small max-temp
```

## Plots

These instructions outline how to generate the plots for:
- Static characterization
- Dynamic characterization

Since the dynamic characterization requires elevated privileges, the analysis needs to run inside the provided docker container.

```
sudo docker pull ghcr.io/binpash/benchmarks:latest
sudo docker run -it --rm -v "$(pwd)":/mnt ghcr.io/binpash/benchmarks:latest bash
```

Then, inside the container run:

```
./infrastructure/generate-plots.sh
```

The plots should be inside the `/mnt` directory in the container, and
the `pwd` directory on the host system. Syntax analysis heatmap will be located at infrastructure/plots/koala-stx-analysis.pdf and dynamic analysis plots will be located in infrastructure/plots/dynamic_analysis.


# Contact

For questions or bug reports, please contact the authors or open an issue on GitHub.
