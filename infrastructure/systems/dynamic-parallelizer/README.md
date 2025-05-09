## hs README

### Overview

`hs` is a system for executing shell scripts out of order. It achieves this by tracing the script's execution, and if an error arises due to speculative execution, the script re-executes the necessary parts to ensure correct outcomes. The project aims to boost the parallel execution of shell scripts, reducing their runtime and enhancing efficiency.

### Security Warning

Being an experimental project, we are currently using sudo to change permission of `/sys/fs/cgroup/cgroup.procs` to 666 (which by default is usually 644).

### Structure

The project's top-level directory contains the following:

- `deps`: Dependencies required by `hs`.
- `docs`: Documentation and architectural diagrams.
- `model-checking`: Tools and utilities for model checking.
- `parallel-orch`: Main orchestration components.
- `pash-spec.sh`: Entry script to initiate the `hs` process.
- `README.md`: This documentation file.
- `report`: Generated reports related to test runs and performance metrics.
- `requirements.txt`: List of Python dependencies.
- `Rikerfile`: Configuration file for Riker.

### Installation

Install `hs` on your Linux-based machine by following these steps:

**Note:** Currently works with `Ubuntu 20.04` or later

1. Navigate to the project directory:
   ```sh
   cd path_to/dynamic-parallelizer
   ```

2. Run the installation script:
   ```sh
   ./scripts/install_deps_ubuntu20.sh
   ```

This script will handle all the necessary installations, including dependencies, try, Riker, and PaSh.

### Running `hs`

The main entry script to initiate `hs` is `pash-spec.sh`. This script sets up the necessary environment and invokes the orchestrator in `parallel-orch/orch.py`. It's designed to accept a variety of arguments to customize its behavior, such as setting debug levels or specifying log files.

Example of running the script:

```bash
./pash-spec.sh [arguments] script_to_speculatively_run.sh
```

**Arguments**:

- `-d, --debug-level`: Set the debugging level. Default is `0`.
- `-f, --log_file`: Define the logging output file. By default, logs are printed to stdout.
- `--sandbox-killing`: Kill any running overlay instances before committing to the lower layer.
- `--env-check-all-nodes-on-wait`: On a wait, check for environment changes between the current node and all other waiting nodes. (not fully functional yet!)

### Testing

To run the provided tests:

```bash
./test/test_orch.sh
```

For in-depth analysis, set the `DEBUG` environment variable to `2` for detailed logs and redirect logs to a file:

```bash
DEBUG=2 ./test/test_orch.sh 2>logs.txt
```

### Contributing and Further Development

Contributions are always welcome! The project roadmap includes extending the architecture to support complete scripts, optimizing the scheduler for better performance, etc.

For a detailed description of possible optimizations, see the [related issues](https://github.com/binpash/dynamic-parallelizer/issues?q=is%3Aopen+is%3Aissue+label%3Aoptimization)

### License

`hs` is licensed under the MIT License. See the `LICENSE` file for more information.
