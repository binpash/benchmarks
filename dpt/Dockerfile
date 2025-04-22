FROM debian:12.7

WORKDIR /benchmarks
COPY ./dpt/ /benchmarks/dpt/
COPY ./infrastructure/ /benchmarks/infrastructure/
COPY ./main.sh /benchmarks/main.sh
COPY ./.git /benchmarks/.git


RUN apt update && apt install -y --no-install-recommends \
    sudo \
    curl \
    wget \
    unzip \
    python3-pip \
    git \
    gpg

RUN useradd -m user && \
    echo "user ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    chown -R user:user /benchmarks

RUN git config --global --add safe.directory /benchmarks

CMD ["bash"]
