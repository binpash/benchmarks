# TERA-Seq Benchmark

> Downloads and sequences genome data

# Building TERA-Seq

Clone and go into this directory.

1. Inputs: `./inputs/inputs.sh`
2. Building image: `docker build -t teraseq-base .` Warning: 177GB Docker image!
3. Running image: `docker run -it --rm teraseq-base /bin/bash`
