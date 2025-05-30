## file-enc

This benchmark compresses and encrypts a collection of input files using gzip and AES-256 encryption.

### Inputs

- `inputs/pcaps/`: A directory containing input files to be processed.

### Running

Two processing steps are applied to each input file:

1. Compression:
   - Each file is compressed using gzip.
   - The result is saved to `outputs/compress_files/<filename>.zip`.

2. Encryption:
   - Each file is encrypted using OpenSSL with AES-256-CBC.
   - The result is saved to `outputs/encrypt_files/<filename>.enc`.

### Validation

Correctness is determined by computing the MD5 hash of each output file and comparing it against a reference hash stored in `hashes/`.

### References

- https://arxiv.org/abs/2012.10206
