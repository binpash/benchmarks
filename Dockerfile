FROM debian:12.7

WORKDIR /benchmarks

COPY . .

RUN apt update && apt install -y --no-install-recommends \
    sudo \
    tcpdump \
    curl \
    wget \
    unzip \
    git \
    samtools \
    minimap2 \
    bcftools \
    python3-pip \
    vim \
    ffmpeg \
    unrtf \
    imagemagick \
    libarchive-tools \
    libncurses5-dev \
    libncursesw5-dev \
    zstd \
    liblzma-dev \
    libbz2-dev \
    zip \
    unzip \
    nodejs \
    tcpdump \
    autoconf \
    automake \
    libtool \
    build-essential \
    gawk \
    pkg-config

RUN pip3 install --break-system-packages \
    scikit-learn \
    kaggle

# Install PaSh
RUN wget https://raw.githubusercontent.com/binpash/pash/main/scripts/up.sh && \
    sh up.sh && \
    export PASH_TOP="$PWD/pash/"

# Test PaSh installation
RUN "$PWD/pash/pa.sh" -c "echo hi"

# Default command
CMD ["bash"]
