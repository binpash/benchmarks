# Overview

The paper makes the following claims on pg. 2 (Comments to AEC reviewers after `:`):

1. **A collection of benchmark programs**: these are POSIX shell programs that perform a variety of tasks.

2. **Permanent and scalable storage of varying inputs**: "permanent" means Zenodo; "scalable" means Brown university servers, which allows for faster download compared to Zenodo. There are three input sizes---`min` take about a few minutes, `small` take a few hours, and `full` take several days---but note _AEC reviewers are encouraged to only try the first two sizes_.

3. **Tooling and automation for reuse**: additional scripts automating dependency setup, input validation, and correctness checking.

4. **Characterization and implications**: additional analysis performing static and dynamic analysis of the benchmark programs.

This artifact targets the following badges:

* [ ] [Artifact available](#artifact-available): Reviewers are expected to confirm that the benchmark programs, their inputs, and automation scripts are all publicly available (about 10 minutes).
* [ ] [Artifact functional](#artifact-functional): Reviewers are expected to confirm sufficient documentation, key components as described in the paper, and execution with `min` inputs (about 20 minutes); execution with `small` inputs (3 hours)is encouraged only after completion of the full artifact evaluation.
* [ ] [Results reproducible](#results-reproducible): Reviewers are expected to reproduce _key_ results of sections 3, 5, 6, and 7 of the paper (2 hours).

# Artifact Available (10 minutes)

The KOALA benchmark suite is:

1. The benchmark code is permanently hosted at: [https://github.com/binpash/benchmarks](https://github.com/binpash/benchmarks)
2. Data are available on _permanent_ (i.e., archival-level durable, but slow) and _scalable_ (i.e., fast, but not archival-level durable) storage:
    1. Permanent arvhival storage on Zenodo, split across multiple DOIs due to Zenodo's max 50GB limit (AEC Reviewers: this is slow, do not try to download---just confirm their existence):
        - [Small-size inputs](https://zenodo.org/records/15361083)
        - [Full inputs (1/5)](https://zenodo.org/records/15367723)
        - [Full inputs (2/5)](https://zenodo.org/records/15368074)
        - [Full inputs (3/5)](https://zenodo.org/records/15368508)
        - [Full inputs (4/5)](https://zenodo.org/records/15368510)
        - [Full inputs (5/5)](https://zenodo.org/records/15368512)
    2. Scalable storage on a Brown University cluster (additional clusters being set up at Stevens and UCLA), accessible via `https://atlas.cs.brown.edu/data`.
3. All scripts required to reproduce the results are available in the [`infrastructure/`](https://github.com/binpash/benchmarks/tree/main/infrastructure) directory of the benchmark repository.

**Please use the Brown links for testing and evaluation, as the Zenodo links take a long time to download.**

# Artifact Functional

## Documentation

[Repo README](https://github.com/binpash/benchmarks)

Per-benchmark documentation is available within each benchmark directory, it
contains information about the computation, input data, and expected output:
[aurpkg](https://github.com/binpash/benchmarks/tree/main/aurpkg),
[bio](https://github.com/binpash/benchmarks/tree/main/bio),
[covid-mts](https://github.com/binpash/benchmarks/tree/main/covid-mts),
[file-enc](https://github.com/binpash/benchmarks/tree/main/file-enc),
[log-analysis](https://github.com/binpash/benchmarks/tree/main/log-analysis),
[makeself](https://github.com/binpash/benchmarks/tree/main/makeself),
[max-temp](https://github.com/binpash/benchmarks/tree/main/max-temp),
[media-conv](https://github.com/binpash/benchmarks/tree/main/media-conv),
[nlp](https://github.com/binpash/benchmarks/tree/main/nlp),
[oneliners](https://github.com/binpash/benchmarks/tree/main/oneliners),
[riker](https://github.com/binpash/benchmarks/tree/main/riker),
[sklearn](https://github.com/binpash/benchmarks/tree/main/sklearn),
[unix50](https://github.com/binpash/benchmarks/tree/main/unix50),
[vps-audit](https://github.com/binpash/benchmarks/tree/main/vps-audit),
[web-index](https://github.com/binpash/benchmarks/tree/main/web-index).

## Completeness

All 15 benchmarks outlined in the paper are available:
[aurpkg](https://github.com/binpash/benchmarks/tree/main/aurpkg),
[bio](https://github.com/binpash/benchmarks/tree/main/bio),
[covid-mts](https://github.com/binpash/benchmarks/tree/main/covid-mts),
[file-enc](https://github.com/binpash/benchmarks/tree/main/file-enc),
[log-analysis](https://github.com/binpash/benchmarks/tree/main/log-analysis),
[makeself](https://github.com/binpash/benchmarks/tree/main/makeself),
[max-temp](https://github.com/binpash/benchmarks/tree/main/max-temp),
[media-conv](https://github.com/binpash/benchmarks/tree/main/media-conv),
[nlp](https://github.com/binpash/benchmarks/tree/main/nlp),
[oneliners](https://github.com/binpash/benchmarks/tree/main/oneliners),
[riker](https://github.com/binpash/benchmarks/tree/main/riker),
[sklearn](https://github.com/binpash/benchmarks/tree/main/sklearn),
[unix50](https://github.com/binpash/benchmarks/tree/main/unix50),
[vps-audit](https://github.com/binpash/benchmarks/tree/main/vps-audit),
[web-index](https://github.com/binpash/benchmarks/tree/main/web-index).

## Exercisability

Follow the instructions below:

**Quickstart:**

To run the benchmarks directly on your machine, without containerization, use the following lines.
Note that the scripts *assume* some basic dependencies that very minimal system setups may not already have (e.g. `git`).

```sh
git clone https://github.com/binpash/benchmarks
cd benchmarks
./kick-tires.sh
```
Expected time: ~2 -- 20 minutes, depending on hardware and network speed.
Runs all benchmarks with minimal inputs (most of the time will be spent on downloading dependencies).

**Benchmark-specific example (nlp):**

```sh
cd benchmarks
./main.sh --min nlp
```

**Configurable system integration:**

To run each of the benchmarks on a system that acts as a shell interpreter, you
can set the `KOALA_SHELL` environment variable.

```sh
KOALA_SHELL="bash --posix" ./main.sh --min max-temp
```

**Bare run:**
To run a benchmark directly on your machine, without containerization, use the following lines.

```
./main.sh --min --bare max-temp
```

**Small inputs:**

To run the benchmarks with small inputs, use the following lines.

```sh
./kick-tires.sh --small
```

# Results Reproducible

The main results presented in the paper are:
1. Results of the static characterization of the benchmarks
2. Results of the dynamic characterization of the benchmarks
3. Results of PCA analysis of (1) the dynamic characterization and (2) the benchmark source code
4. Results of the application of the benchmark suite to a set of prior shell optimization tools

## Static Characterization

To generate the static characterization, run the following commands:

```sh
# TODO
```


```
mkdir -p /tmp/plots
sudo docker pull ghcr.io/binpash/benchmarks:latest
sudo docker run -it --rm -v "/tmp/plots":/mnt ghcr.io/binpash/benchmarks:latest bash
```

Then, inside the container run:

```
./infrastructure/generate-plots.sh
```

The plots should be inside the `/mnt` directory in the container, and
the `/tmp/plots` directory on the host system.
Syntax analysis heatmap will be located at `infrastructure/plots/koala-stx-analysis.pdf` and dynamic analysis plots will be located in `infrastructure/plots/dynamic_analysis`.

## Dynamic Characterization

To generate the dynamic characterization, run the following commands:

```sh
# TODO
```

## PCA Analysis

The second part of the PCA analysis involves sending the source code of the benchmarks to a remote embedding model using OpenAI's API.
For convenience and cost concerns, we provide the results the embedding model in the `infrastructure/plots/pca` directory.

```sh
# TODO
```

# Contact

For questions or bug reports, please contact the authors or open an issue on GitHub.
