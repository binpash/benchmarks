## README for `hs` Test Suite

### Overview
This directory contains the test suite for `hs`. The main test script is `test_orch.sh`, which automates the process of running various tests on the `hs` and `bash` to ensure consistency and correctness.

### Directory Structure

- **test_scripts**: Contains the individual test scripts.
- **misc**: Contains utility scripts used by the test cases.
- **output_bash**: Directory to save the output of scripts executed by `bash`.
- **output_orch**: Directory to save the output of scripts executed by `hs`.
- **results**: Stores the result status and logs for each test.
- **parse_cmd_repetitions.py**: Python script to parse command repetitions from the `hs` logs.

### Main Test Script (`test_orch.sh`)

The main test script `test_orch.sh` starts by setting up environment variables and directories. It then proceeds to define utility functions:

- `cleanup()`: Removes cache and clears output directories.
- `test_repetitions()`: Validates the repetition of commands using `parse_cmd_repetitions.py`.
- `run_test()`: Executes a given test for both `bash` and `hs` and compares the outputs.
- Various test functions, e.g., `test_single_command()`, `test_local_vars_1()`, etc.

Finally, it runs the set of defined tests, provides a summary of the results, and outputs logs for both passed and failed tests.

### Running Tests

To run all tests:
```
./test_orch.sh
```

To run specific tests:
```
./test_orch.sh [testname]
```

Before running your scripts, you can set the DEBUG environment variable to provide detailed logging information. Assign a value of 2 to DEBUG to get the most detailed logs.

```bash
export DEBUG=2
```

Since the logs are printed to stderr, you can redirect them to a file to facilitate easier analysis:

```bash
./test_orch.sh [test_name] 2>logs.txt
```

### Test Results

At the end of execution, a summary is presented:

1. List of tests that produced identical outputs on both `bash` and `hs`.
2. List of tests that produced non-identical outputs.
3. Overall summary indicating the number of tests passed.

The detailed logs for passed and failed tests can be found in the `results` directory.

### Adding More Tests

If you would like to expand the test suite by adding more tests, follow these guidelines:

1. **Create Test Script**: Write a new Bash script that performs the desired test. For example, if you wish to test a new functionality named `test_XXX.sh`, create a file with that name under the `test_scripts` directory. <br><br>
Utilize utility scripts from the `misc` directory, like `sleep_and_grep.sh`, to help maintain a modular design. This way, changes made to utility functions can propagate across multiple tests.

2. **Update Main Test Suite**:
    In the main test suite (shown above), add a new function named similarly to your test script. This function should prepare any required input files and run your test script. For example:

    ```bash
    test_XXX()
    {
        local shell=$1
        # Setup input data here if required
        $shell $TEST_SCRIPT_DIR/test_XXX.sh
    }
    ```

3. **Integrate into Test Runner**:
    Add a call to your `run_test` function with your new function as the argument. This ensures it's part of the suite when no specific test names are given. Add this before `if [ "$#" -eq 0 ]; then`, for instance:

    ```bash
    run_test test_XXX
    ```

    Optionally, you can provide expected repetition values as a second argument if required by the test.
