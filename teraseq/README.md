# TeraSeq Benchmark

> Downloads and sequences genome data

# Running TeraSeq

Manual:

1. Clone this repo into a disk that has space. Also ensure that the Docker data root is on a disk that has space.
2. `cd benchmarks/biomed/teraseq`
3. `docker build -t teraseq-bench .` Warning: 127GB Docker image!
4. `cd 5TERA` (Replace `5TERA` with a directory of your choosing)
5. Download: `docker run --rm -v .:/root/TERA-Seq_manuscript/samples teraseq-bench /bin/bash /root/TERA-Seq_manucript/samples/download.sh` (Replace `5TERA` with a directory of your choosing)
6. Run: `time docker run --rm -v .:/root/TERA-Seq_manuscript/samples teraseq-bench /bin/bash /root/TERA-Seq_manucript/samples/run.sh`
7. Compute hash: (bash-hashes.log)
8. Change command in step 6 accordingly to use hs. `...`
9. Compute hash: (bash-hashes.log)
10. Verify: `./verification.sh`
11. Clean up all files (including downloads): `./cleanup.sh -a`
