## makeself

This benchmark runs a suite of shell script tests to evaluate the behavior and correctness of self-extracting archive operations.

### Inputs

- This benchmark has no external inputs.

### Running

1. Each test script under `makeself/test/` is executed using the shell under test.
2. Test output is written to `test_results.log` files inside the corresponding test directories.
3. A summary of test pass/fail status is recorded in `run_results.log`.

### Validation

Correctness is determined by checking for the presence of any failed tests in `run_results.log`.

### References

- https://makeself.io
