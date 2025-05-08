## nlp

This benchmark performs a suite of text-processing analyses over a corpus of digitized books from Project Gutenberg.

### Inputs

- `inputs/pg/`: A directory of plain-text book files.

### Running

Each script in the `scripts/` directory performs a distinct natural language processing task on a subset of the book files. Tasks include:

- Token-level and type-level word counting
- Frequency analysis of bigrams, trigrams, and syllables
- Pattern detection (e.g., uppercase words, words without vowels, common prefixes/suffixes)
- Text comparisons between specific documents (e.g., Genesis and Exodus)

The output of each script is written to a separate subdirectory under `outputs/`, one file per input document.

### Validation

Correctness is determined by computing the SHA-256 hash of each output file and comparing it against a reference hash stored in `hashes/`.

### References

- https://web.stanford.edu/class/cs124/kwc-unix-for-poets.pdf
