#!/bin/bash

if [[ -z "${1:-}" ]]; then
    echo "Usage: $0 <benchmark-name> [flags]"
    exit 1
fi
BENCHMARK=$(basename "$1")
shift
export BENCHMARK

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_SCRIPT

if [[ -z "$BENCHMARK" ]]; then
    echo "Usage: $0 <benchmark-name> [flags]"
    exit 1
fi

BENCHMARK_DIR="./$BENCHMARK"
DOCKERFILE="$BENCHMARK_DIR/Dockerfile"
IMAGE_NAME="koala-$BENCHMARK"

if [[ ! -d "$BENCHMARK_DIR" ]]; then
    echo "Error: Benchmark directory '$BENCHMARK_DIR' not found."
    exit 1
fi

echo "Generating Dockerfile for benchmark: $BENCHMARK"

# Create a new Dockerfile
cat <<EOF >"$DOCKERFILE"
FROM debian:12.7

WORKDIR /benchmarks
COPY ./$BENCHMARK/ /benchmarks/$BENCHMARK/
COPY ./infrastructure/ /benchmarks/infrastructure/
COPY ./main.sh /benchmarks/main.sh
COPY ./.git /benchmarks/.git


RUN apt update && apt install -y --no-install-recommends \\
    sudo \\
    curl \\
    wget \\
    unzip \\
    python3-pip \\
    git \\
    gpg

RUN useradd -m user && \\
    echo "user ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \\
    chown -R user:user /benchmarks

RUN git config --global --add safe.directory /benchmarks

USER user

CMD ["bash"]
EOF
#TODO fix ownership issues ^^

echo "Building Docker image: $IMAGE_NAME..."
docker build -t "$IMAGE_NAME" -f "$DOCKERFILE" . || {
    echo "Docker build failed"
    exit 1
}

echo "Running Docker image: $IMAGE_NAME..."
docker run --cap-add NET_ADMIN --cap-add NET_RAW --rm -it "$IMAGE_NAME"

echo "Docker image $IMAGE_NAME built and run successfully."
echo "Cleaning up..."
docker rmi "$IMAGE_NAME" || {
    echo "Failed to remove Docker image $IMAGE_NAME"
    exit 1
}
rm -f "$DOCKERFILE" || {
    echo "Failed to remove Dockerfile $DOCKERFILE"
    exit 1
}
echo "Docker image $IMAGE_NAME and $DOCKERFILE removed successfully."
