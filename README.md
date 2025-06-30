# The Koala Benchmark Suite

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

## Quick Setup

The suite can be obtained using the following ways:

* Run `curl up.kben.sh | sh` from your terminal, or
* [Clone the repo](https://github.com/kbensh/koala) and run `cd koala && ./setup.sh`, or
* Fetch the Docker container by running `docker pull ghcr.io/kbensh/koala:latest`, or
* Build the Docker [container from scratch](#environment--setup-notes).

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
