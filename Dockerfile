FROM ubuntu:20.04

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies that are always required
RUN apt-get update && apt-get install -y \
    bash \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory to the root folder containing all benchmarks
WORKDIR /benchmarks

COPY . /benchmarks
RUN chmod +x /benchmarks/*/*.sh

CMD cd "$BENCHMARK_FOLDER" && \
    ./dependencies && \
    ./inputs.sh && \
    ./run.sh && ./verify.sh && ./cleanup.sh || ./cleanup.sh
