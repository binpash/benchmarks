FROM debian:12.7

WORKDIR /benchmarks

RUN apt update && apt install -y --no-install-recommends \
    sudo \
    curl \
    wget \
    unzip \
    python3-pip \
    git \
    gpg

RUN useradd -m user && \
    echo "user ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
COPY . .

RUN chown -R user:user /benchmarks
RUN chmod +x /benchmarks/main.sh
USER user

# Run install.sh for each benchmark
RUN set -eux; \
    for bench in /benchmarks/*; do \
        if [ -f "$bench/install.sh" ]; then \
            echo "Running install.sh in $bench"; \
            chmod +x "$bench/install.sh"; \
            bash "$bench/install.sh"; \
        fi; \
    done

RUN git config --global --add safe.directory /benchmarks

CMD ["bash"]
