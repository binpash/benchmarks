#!/bin/bash

export BENCHMARK=$(basename "$1")
shift

if [[ -z "$BENCHMARK" ]]; then
    echo "Usage: $0 <benchmark-name> [flags]"
    exit 1
fi

BENCHMARK_DIR="./$BENCHMARK"
DOCKERFILE="$BENCHMARK_DIR/Dockerfile"

if [[ ! -d "$BENCHMARK_DIR" ]]; then
    echo "Error: Benchmark directory '$BENCHMARK_DIR' not found."
    exit 1
fi

echo "Generating Dockerfile for benchmark: $BENCHMARK"

# Create a new Dockerfile
cat <<EOF > "$DOCKERFILE"
FROM debian:12.7

WORKDIR /benchmarks
COPY ./$BENCHMARK /benchmarks/$BENCHMARK
COPY ./.git /benchmarks/.git
COPY ./main.sh /benchmarks/main.sh
COPY ./infrastructure /benchmarks/infrastructure


RUN apt update && apt install -y --no-install-recommends \\
    sudo \\
    tcpdump \\
    curl \\
    wget \\
    unzip \\
    bcftools \\
    python3-pip \\
    vim \\
    git \\
    gpg

RUN useradd -m user && \\
    echo "user ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \\
    chown -R user:user /benchmarks

RUN git config --global --add safe.directory /benchmarks

CMD ["bash"]

EOF

IMAGE_NAME="koala-$BENCHMARK"
echo "Building Docker image: $IMAGE_NAME..."
docker build -t "$IMAGE_NAME" -f "$DOCKERFILE" . || { echo "Docker build failed"; exit 1; }
echo "Running Docker image: $IMAGE_NAME..."
docker run --cap-add NET_ADMIN --cap-add NET_RAW --rm -it $IMAGE_NAME
