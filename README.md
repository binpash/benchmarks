# The Koala Benchmark Suite
[Benchmarks](#benchmarks) | [Quick Setup](#quick-setup) | [More Info](#more-info) | [License](#license)

> If you find any *issues*, [open a GitHub issue](https://github.com/binpash/benchmarks/issues/new/choose).

The Koala benchmark suite enables characterization of performance-oriented research targeting the POSIX shell. It consists of 14 real-world program sets from diverse domains ranging from CI/CD and AI/ML to biology and the humanities. They are accompanied by real inputs that facilitate small- and large-scale performance characterization and varying opportunities for optimization.

If any aspect of Koala is useful, please cite the [ATC'25 Koala paper](https://www.usenix.org/conference/atc25/presentation/lamprou):
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

As part of the [ATC'25 Artifact Evaluation process](https://www.usenix.org/conference/atc25/call-for-artifacts), the Koala frozen [`atc25-ae branch`](https://github.com/kbensh/koala/tree/atc25-ae) received all three badges—artifact *Available*, *Functional*, and *Reproduced*.

## Benchmarks

| Benchmark   | Description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| `analytics` | processes real-world network logs to extract and summarize key events.          |
| `bio`       | performs genomic and transcriptomic analysis using population and RNA-seq data. |
| `ci-cd`     | builds and tests open-source software projects.                             |
| `covid`     | analyzes public transit activity during the covid-19 pandemic.              |
| `file-mod`  | compresses, encrypts, and converts various file formats.                    |
| `inference` | runs media-related inference tasks using large foundation models.           |
| `ml`        | implements a full machine learning pipeline using scikit-learn.             |
| `nlp`       | processes books using shell-based nlp pipelines from unix for poets.        |
| `oneliners` | executes classic and modern one-liner shell pipelines.                      |
| `pkg`       | builds aur packages and analyzes npm packages for permissions.              |
| `repl`      | performs security auditing and git-based development workflow replay.       |
| `unixfun`   | solves unix text-processing problems from the 50-year anniversary challenge.|
| `weather`   | computes and visualizes historical weather statistics.                      |
| `web-search`| implements crawling, indexing, and querying of wikipedia data.              |

## Quick Setup

Koala can be obtained using the following ways:

* Run `curl up.kben.sh | sh` from your terminal,
* [Clone the repo](https://github.com/kbensh/koala) and run `cd koala && ./setup.sh`,
* Fetch the oficial Docker container by running `docker pull ghcr.io/kbensh/koala:latest`, or
* Build the Docker [container from scratch](#environment--setup-notes).

## More info

**Setup and configuration:** Extended instructions on Koala's setup and configuration can be found in the [INSTRUCTIONS](./INSTRUCTIONS.md) file.

**Inputs and dependencies:** Information on Koala's multiple inputs (e.g., min, small, and full)  as well as the dependencies of all benchmarks, can be found in the [INSTRUCTIONS](./INSTRUCTIONS.md) file.

**Contributions, always welcomed!** Further details on how to contribute to the Koala benchmark project, take a look at the [CONTRIBUTING](./CONTRIBUTING.md) file.

# License

The Koala Benchmarks are licensed under the MIT License. See the LICENSE file for more information.
