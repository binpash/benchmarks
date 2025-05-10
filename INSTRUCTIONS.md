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

Confirm that the benchmark programs, their inputs, and automation scripts are all publicly available:

1. The benchmark code is permanently hosted at: [https://github.com/binpash/benchmarks](https://github.com/binpash/benchmarks)

3. Data are available on _permanent_ (i.e., archival-level durable, but slow) and _scalable_ (i.e., fast, but not archival-level durable) storage:

    * Permanent archival storage on Zenodo, split across multiple DOIs due to Zenodo's max 50GB limit (AEC Reviewers: _this is slow, do not try to download_—just confirm their existence): [code & software dependencies](https://zenodo.org/records/15377017), [`small` inputs](https://zenodo.org/records/15361083); `full` in 5 parts:
   [1](https://zenodo.org/records/15367723),
   [2](https://zenodo.org/records/15368074)
   [3](https://zenodo.org/records/15368508),
   [4](https://zenodo.org/records/15368510),
   [5](https://zenodo.org/records/15368512).

* Scalable storage on a Brown University cluster (additional clusters being set up at Stevens and UCLA), accessible via `https://atlas.cs.brown.edu/data`.

3. Additional scripts are available in the [`infrastructure/`](https://github.com/binpash/benchmarks/tree/main/infrastructure) directory of the repository.

> AEC Reviewers: From this point on, scripts use the Brown links, as Zenodo is significantly slower.

# Artifact Functional (20 minutes, optionally 3 hours)

Confirm sufficient documentation, key components as described in the paper, and execution with min inputs (about 20 minutes):

* Documentation: [Repo README](https://github.com/binpash/benchmarks) and benchmark-specific documentation available within each benchmark directory,
contains information about each computation, input data, and expected output:
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

* Key components: 15 benchmarks (i.e.,
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
[web-index](https://github.com/binpash/benchmarks/tree/main/web-index)) and several scripts in [`infrastructure/`](https://github.com/binpash/benchmarks/tree/main/infrastructure).

* Exercisability: Instructions below set up an Debian-based container and run _all_ benchmarks on `min` inputs (`all.sh`) or run specific benchmarks.

> At this point, **run `git clone https://github.com/binpash/benchmarks` and `cd benchmarks`.**

**Quickstart: Running a single benchmark (e.g., `nlp`):** To quickly execute a specific benchmark such as `nlp`, invoke the top-level `main.sh` script—which will set up a Debian-based container image, install dependencies, download benchmark-specific `min` inputs, and execute the benchmark (2 minutes):

```sh
./main.sh nlp
```

The reason this is easier to evaluate in a container is that some scripts *assume* that some dependencies on the *host system* (e.g. `git`). The scripts also assume a docker executable in the PATH, which can be configured with the `KOALA_CONTAINER_CMD`.

**Complete run: Running all benchmarks:** To execute the full set of benchmarks, invoke the top-level `all.sh` script—which essentially runs the previous command in a loop (20 minutes):

```sh
cd benchmarks
./all.sh
```

> **Not recommended on the first run:**
> 
> Both scripts above take additional optional parameters: Using `--small` will execute the benchmark(s) with the `small` inputs (about 3 hours) and `--bare` will execute the benchmarks on the host system (which might not satisfy dependencies). To replace the shell interpreter, `set` the `KOALA_SHELL` variable — e.g., `KOALA_SHELL="bash --posix"`.


# Results Reproducible

The key results of the analysis are the following:

1. Principal Component Analysis (§3, Fig. 2);
2. Syntactic Characterization & Analysis (§5, Fig. 4 & 5);
3. Dynamic Characterization & Analysis (§6, Fig. 6)

Results 1–3 are easy to reproduce, by running specific scripts (a part of PCA depends on embeddings, which incurs some financial costs due to OpenAI); PCA depends on the results of the analyses, so we start with no. 2–3 before returning to no. 1.

**Preparation:** The dynamic analysis requires running Koala in a pre-set container due to elevated privileges:

```sh
mkdir -p /tmp/plots
sudo docker pull ghcr.io/binpash/benchmarks:latest
sudo docker run -it --rm -v "/tmp/plots":/tmp/plots ghcr.io/binpash/benchmarks:latest bash
./setup.sh # Inside the container
```

**Static Characterization:** To generate the static characterization, run the following commands:

```sh
./infrastructure/static-analysis.sh /tmp/plots
```

This will place the static analysis heat-map in the `/tmp/plots` directory on the host system.

**Dynamic Characterization:** To generate the dynamic characterization, run the following commands:
*Note*: This step will run the entire suite again, now having tracing enabled (~3 hours).

```sh
./infrastructure/dynamic-analysis.sh --small
```

**PCA Analysis:** The second part of the PCA analysis involves sending the source code of the benchmarks to a remote embedding model using OpenAI's API.
For convenience and cost concerns, we provide the results the embedding model in `infrastructure/data/embeddings.csv` directory.
This file can also be generated by running `python3 infrastructure/do_embedding.py`.

```sh
./infrastructure/pc-analysis.sh /tmp/plots/pca.pdf
```

## Optional: Applying benchmarks to various systems (1–3 days)

This section is about setting up and running other systems on the Koala benchmarks. Crucially:

* The difficulty of evaluating prototype systems using Koala depends on these systems (i.e., not our contribution) as well as access to specific hardware (e.g., large multiprocessors).
* The authors of these systems are free to make any modifications to the benchmarks, the inputs, or the scripts according to the contributions claimed by _their_ systems, not Koala

Therefore, this evaluation outside the scope of the Koala artifact evaluation.

# Contact

For questions or bug reports, please contact the authors or open an issue on GitHub.
