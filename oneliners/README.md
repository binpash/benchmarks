## oneliners

This benchmark runs a collection of small, single-purpose shell scripts over large text files to evaluate data processing tasks such as sorting, frequency counting, and regex-based pattern matching.

### Inputs

- `inputs/`: A collection of text files and structured data files (`logs-popcount-org.txt`, `dict.txt`).

### Running

Each script processes one input file using a Unix pipeline. Examples include:

- `sort.sh`: Sorts all lines of the input.
- `top-n.sh`: Computes the top 1000 most frequent words.
- `wf.sh`: Outputs word frequencies sorted by count.
- `diff.sh`: Compares case-normalized versions of a stream.
- `nfa-regex.sh`: Applies a complex regular expression over the input.
- `uniq-ips.sh`: Extracts unique IP addresses from logs.

Outputs are stored as `outputs/<script>.out`.

### Validation

Correctness is determined by computing the SHA-256 hash of each output file and comparing it against a reference hash stored in `hashes/`.
