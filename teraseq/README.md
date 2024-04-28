# TERA-Seq Benchmark

> Downloads and sequences genome data

# Running TERA-Seq

1. Clone this repo into a disk that has space 600GB of free space should be good enough for the full TERA-Seq run without intermediate cleanups. Inputs total around 46GB. Also ensure that the Docker data root is on a disk that has space.
2. `cd benchmarks/teraseq`
3. `docker build -t teraseq-data .` Warning: 100GB+ Docker image!
4. `docker run -it --rm teraseq20-data /bin/bash`

# Building Tools on top of TERA-Seq

TBD
