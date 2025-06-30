## ci-cd

This benchmark suite builds and tests multiple real-world software systems from
source (`vim`, `redis`, `sqlite`, `lua`, `memcached`, `xz`) and evaluates
the behavior and correctness of self-extracting archive utility `makeself`.

### Inputs

- Software Builds:
  - Each subbenchmark clones a pinned version of the corresponding repository into:
    ```
    inputs/scripts/<project>/dev/
    ```

- Makeself:
  - No external inputs are required. Test scripts are included under:
    ```
    makeself/test/
    ```

### Running

1. Build and Test Software:
   - For each software project:
     - Clone the repository.
     - Execute the hand-crafted, dependency-aware build script to simulate a complete build pipeline.
     - Run the associated `validate.sh` script to ensure correct build and behavior.

2. Validate Makeself Archives:
   - For each makeself test:
     - Execute the test script under the shell under test.
     - Write detailed output logs to:
       ```
       makeself/test/<test_name>/test_results.log
       ```
     - Collect a summary of test pass/fail status in:
       ```
       makeself/run_results.log
       ```

### Validation

Validation is performed in two parts:

- Software Builds:
  - Each `validate.sh` script checks that:
    - A specific binary was successfully built.
    - The binary runs and produces the expected output.
  - Example checks:
    - `vim` opens a file and writes output.
    - `redis-cli` returns a version string.

- Makeself Archives:
  - The benchmark verifies that all self-extracting archive tests pass by confirming no failures are recorded in `makeself/run_results.log`.

### References

- [makeself.io](https://makeself.io)
- [riker](https://github.com/curtsinger-lab/riker)
