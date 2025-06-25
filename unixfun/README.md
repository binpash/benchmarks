## unix50

This benchmark replicates the classic "Unix for Beginners" exercises and
trivia-style queries, designed to test proficiency in Unix shell utilities by
operating over structured historical, biographical, and technical text
datasets.

### Inputs

- `inputs/`: A collection of `.txt` files organized by chapter or question number.

### Running

Each script corresponds to a numbered task and operates over its associated input file. Scripts perform classic Unix data wrangling, such as:

- Extracting fields via `cut` or `awk`
- Word frequency analysis via `sort`, `uniq`, and `wc`
- Pattern extraction with `grep`, `tr`, and regular expressions

Outputs are stored in `outputs/<script_number>.out`.

### Validation

SHA-256 hashes of outputs are compared against known-good values in the `hashes/` directory for each input size (`min`, `small`, or `full`) to ensure correctness.

### References

- https://unix50.org/
