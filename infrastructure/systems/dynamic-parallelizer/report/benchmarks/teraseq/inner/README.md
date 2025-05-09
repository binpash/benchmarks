# TeraSeq Benchmark

> Downloads and sequences genome data

# Running TeraSeq

Manual:

1. Clone this repo into a disk that has space. Also ensure that the Docker data root is on a disk that has space.
2. `cd dynamic-parallelizer/report/benchmarks/teraseq`
3. `docker build -t teraseq-data .` Warning: 128GB Docker image!
4. `cd 5TERA` (Replace `5TERA` with a directory of your choosing)
5. Download: `docker run --rm -v .:/root/TERA-Seq_manuscript/samples teraseq-data /bin/bash /root/TERA-Seq_manucript/samples/download.sh` (Replace `5TERA` with a directory of your choosing)
6. Run: `time docker run --rm -v .:/root/TERA-Seq_manuscript/samples teraseq-data /bin/bash /root/TERA-Seq_manucript/samples/run.sh`
7. Compute hash: (bash-hashes.log)
8. Change command in step 6 accordingly to use hs. `...`
9. Compute hash: (bash-hashes.log)
10. Verify: `./verification.sh`
11. Clean up all files (including downloads): `./cleanup.sh -a`

# Building Docker Container with hs
1. `cd dynamic-parallelizer/report/benchmarks/teraseq`
2. `docker build -t teraseq-data .`
3. `cd ../../..`
4. `docker build -t teraseq-hs . -f report/benchmarks/teraseq/Dockerfile.hs`

We build the hs image on top of the data image to allow for quick iteration on our tool.
