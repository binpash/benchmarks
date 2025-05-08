## riker

This benchmark suite builds and tests multiple real-world software systems from source (`vim`, `redis`, `sqlite`, `lua`, `memcached`, `xz`) using hand-crafted, dependency-aware scripts that simulate a full build pipeline.

### Inputs

Each subbenchmark clones a pinned version of the corresponding repository into `input/scripts/<project>/dev/`.

### Validation

Each benchmark validates its results using an associated `validate.sh` script. These usually check that:
	- A specific binary was successfully built.
	- The binary runs and produces the expected output.

Example checks:
	- `vim` opens a file and writes output.
	- `redis-cli` returns a version string.
