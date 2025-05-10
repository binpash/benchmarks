# Overview

The paper makes the following claims on pg. 2 (Comments to AEC reviewers after `:`):

1. **A collection of benchmark programs**: these are POSIX shell programs that perform a variety of tasks.

2. **Permanent and scalable storage of varying inputs**: "permanent" means Zenodo; "scalable" means Brown university servers, which allows for faster download compared to Zenodo. There are three input sizes---`min` take about a few minutes, `small` take a few hours, and `full` take several days---but note _AEC reviewers are encouraged to only try the first two sizes_.

3. **Tooling and automation for reuse**: additional scripts automating dependency setup, input validation, and correctness checking.

4. **Characterization and implications**: additional infrastructure performing static and dynamic program analysis of the benchmark programs.

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

* Scalable storage on a Brown University cluster (additional clusters being set up at Stevens and UCLA), accessible via `https://atlas.cs.brown.edu/data` (see [the full table for _all_ inputs in **Appendix I: All Inputs** below](https://github.com/binpash/benchmarks/blob/main/INSTRUCTIONS.md#appendix-i-all-inputs).)

3. Additional scripts are available in the [`infrastructure/`](https://github.com/binpash/benchmarks/tree/main/infrastructure) directory of the repository.

> AEC Reviewers: From this point on, scripts use the Brown links, as Zenodo is significantly slower.

# Artifact Functional (20 minutes, optionally 3 hours)

Confirm sufficient documentation, key components as described in the paper, and execution with min inputs (about 20 minutes):

* Documentation: The top-level [README](https://github.com/binpash/benchmarks) file and benchmark-specific documentation available within each benchmark directory,
containing information about each computation, input data, and expected output:
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
./all.sh
```

> **Not recommended on the first run:**
> 
> Both scripts above take additional optional parameters: Using `--small` will execute the benchmark(s) with the `small` inputs (about 3 hours) and `--bare` will execute the benchmarks on the host system (which might not satisfy dependencies). To replace the shell interpreter, `set` the `KOALA_SHELL` variable — e.g., `KOALA_SHELL="bash --posix"`.


# Results Reproducible (about 3 hours)

The key results of the analysis are the following:

1. Principal Component Analysis (§3, Fig. 2);
2. Syntactic Characterization & Analysis (§5, Fig. 4 & 5);
3. Dynamic Characterization & Analysis (§6, Fig. 6).

Most of these results are easily reproducible by running specific scripts in the repo. A part of PCA depends on embeddings, whose calculation incurs some financial costs due to OpenAI—thuse we have cashed these embeddings in the repo. Although summarized earlier in the paper (§3), PCA depends on the results of the static and dynamic analyses, so we start with no. 2 and 3 before returning to no. 1.

**Preparation:** The dynamic analysis requires running Koala in a pre-set container due to elevated privileges:

```sh
mkdir -p /tmp/plots
sudo docker pull ghcr.io/binpash/benchmarks:latest
sudo docker run -it --rm -v "/tmp/plots":/tmp/plots ghcr.io/binpash/benchmarks:latest bash
./setup.sh # Inside the container
```

**Static characterization:** To generate the static characterization, run the following commands:

```sh
./infrastructure/static-analysis.sh /tmp/plots
```

This will place the heatmap plot showing the results of the static analysis in the `/tmp/plots` directory on the host system.

**Dynamic characterization:** To generate the dynamic characterization, run the following commands:
*Note*: This step will run the entire suite again, now having tracing enabled (~3 hours).

```sh
./infrastructure/dynamic-analysis.sh --small
```

**Principal component analysis:** Part of the PCA analysis involves sending the source code of the benchmarks to a remote embedding model using OpenAI's API; for convenience and cost concerns, we provide the results the embedding model in `infrastructure/data/embeddings.csv` (generatable by running `python3 infrastructure/do_embedding.py`).

```sh
./infrastructure/pc-analysis.sh /tmp/plots/pca.pdf
```

# Optional: Applying benchmarks to various systems (1–3 days)

This section is about setting up and running other systems on the Koala benchmarks. Crucially:

* The difficulty of evaluating prototype systems using Koala depends on these systems (i.e., not our contribution) as well as access to specific hardware (e.g., large multiprocessors).
* The authors of these systems are free to make any modifications to the benchmarks, the inputs, or the scripts according to the contributions claimed by _their_ systems, not Koala.

Therefore, this evaluation outside the scope of the Koala artifact evaluation.

# Contact

For questions or bug reports, please contact `koala@brown.edu` or open an issue on GitHub.

# Appendix I: All Inputs

The table below contains all links to the inputs. AEC reviewers: Some of these inputs are extremely large and hosted on low-bandwidth permanent-storage services such as Zenodo. To confirm existence, simply click on the URL to to start a download and then stop downloading a few seconds later.

| Benchmark     | Min                                                                                                 | Small                                                                                                            | Full                                                                                                                                  |
|---------------|-----------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| `aurpkg`        | [Github Repo](https://atlas.cs.brown.edu/data/packages_min)                                         | [Brown](https://atlas.cs.brown.edu/data/packages) [Zenodo](https://zenodo.org/records/15361083)                  | [Brown](https://atlas.cs.brown.edu/data/packages) [Zenodo](https://zenodo.org/records/15367723)                                       |
| `bio`           | [Github Repo](https://github.com/binpash/benchmarks/tree/main/bio/min_inputs)                       | [Brown](https://atlas.cs.brown.edu/data/bio/small) [Zenodo](https://zenodo.org/records/15361083)                 | [Brown](https://atlas.cs.brown.edu/data/bio/full) [Zenodo](https://zenodo.org/records/15367723)                                       |
| `covid-mts`     | [Github Repo](https://github.com/binpash/benchmarks/tree/main/covid-mts/min_inputs)                 | [Brown](https://atlas.cs.brown.edu/data/covid-mts/in_small.csv.gz) [Zenodo](https://zenodo.org/records/15361083) | [Brown](https://atlas.cs.brown.edu/data/covid-mts/in_full.csv.gz) [Zenodo](https://zenodo.org/records/15368074)                       |
| `file-enc`      | [Github Repo](https://github.com/binpash/benchmarks/tree/main/file-enc/min_inputs)                  | [Brown](https://atlas.cs.brown.edu/data/pcaps.zip) [Zenodo](https://zenodo.org/records/15361083)                 | [Brown](https://atlas.cs.brown.edu/data/pcaps_large.zip) [Zenodo](https://zenodo.org/records/15368510)                                |
| `log-analysis`  | [Github Repo](https://github.com/binpash/benchmarks/tree/main/log-analysis/min_inputs)              | [Brown](https://atlas.cs.brown.edu/data/pcaps.zip) [Zenodo](https://zenodo.org/records/15361083)                 | [Brown](https://atlas.cs.brown.edu/data/pcaps_large.zip) [Zenodo](https://zenodo.org/records/15368510)                                |
| `makeself`      | No inputs                                                                                           | No inputs                                                                                                        | No inputs                                                                                                                             |
| `max-temp`      | [Github Repo](https://github.com/binpash/benchmarks/tree/main/max-temp/min_inputs)                  | [Brown](https://atlas.cs.brown.edu/data/max-temp/noaa/) [Zenodo](https://zenodo.org/records/15361083)            | [Brown](https://atlas.cs.brown.edu/data/max-temp/noaa/) [Zenodo](https://zenodo.org/records/15368510)                                 |
| `media-conv`    | [Github Repo](https://github.com/binpash/benchmarks/tree/main/media-conv/min_inputs/jpg_min/jpg)    | [Brown](https://atlas.cs.brown.edu/data/media-conv/inputs) [Zenodo](https://zenodo.org/records/15361083)         | [Brown](https://atlas.cs.brown.edu/data/media-conv/inputs) [Zenodo](https://zenodo.org/records/15368510)                              |
| `nlp`           | [Brown](https://atlas.cs.brown.edu/data/gutenberg)                                                  | [Brown](https://atlas.cs.brown.edu/data/gutenberg/) [Zenodo](https://zenodo.org/records/15361083)                | [Brown](https://atlas.cs.brown.edu/data/gutenberg/) [Zenodo](https://zenodo.org/records/15368510)                                     |
| `oneliners`     | [Brown](https://atlas.cs.brown.edu/data/dummy/)                                                     | [Brown](https://atlas.cs.brown.edu/data/dummy/) [Zenodo](https://zenodo.org/records/15361083)                    | [Brown](https://atlas.cs.brown.edu/data/dummy/) [Zenodo](https://zenodo.org/records/15368512)                                         |
| `riker`         | [Github Repo](https://atlas.cs.brown.edu/data/riker)                                                | [Brown](https://atlas.cs.brown.edu/data/riker) [Zenodo](https://zenodo.org/records/15361083)                     | [Brown](https://atlas.cs.brown.edu/data/riker) [Zenodo](https://zenodo.org/records/15368512)                                          |
| `sklearn`       | [Github Repo](https://github.com/binpash/benchmarks/tree/main/sklearn/inputs/covertype)             | [Brown](https://atlas.cs.brown.edu/data/sklearn/) [Zenodo](https://zenodo.org/records/15361083)                  | [Brown](https://atlas.cs.brown.edu/data/sklearn/) [Zenodo](https://zenodo.org/records/15368512)                                       |
| `unix50`        | [Github Repo](https://atlas.cs.brown.edu/data/unix50)                                               | [Brown](https://atlas.cs.brown.edu/data/unix50/small) [Zenodo](https://zenodo.org/records/15361083)              | [Brown](https://atlas.cs.brown.edu/data/unix50/large) [Zenodo](https://zenodo.org/records/15368512)                                   |
| `vps-audit`     | No inputs                                                                                           | No inputs                                                                                                        | No inputs                                                                                                                             |
| `web-index`     | [Brown](https://atlas.cs.brown.edu/data/wikipedia_min.tar.gz)                                       | [Brown](https://atlas.cs.brown.edu/data/wikipedia_small.tar.gz) [Zenodo](https://zenodo.org/records/15361083)    | [Brown](https://atlas.cs.brown.edu/data/wikipedia.tar.gz) [Zenodo](https://zenodo.org/records/15368512)                               |


