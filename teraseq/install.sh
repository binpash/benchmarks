#!/usr/bin/env bash
set -euo pipefail

TOP=$(git rev-parse --show-toplevel)
benchmark_dir="${TOP}/teraseq"

# install.sh: Installs system-wide dependencies for the TERA-Seq pipeline

# 1. Install OS packages
apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    python3 \
    python3-dev \
    python3.11-dev \
    python3-all-dev \
    python3-pip \
    perl \
    cpanminus \
    libdbi-perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    openjdk-17-jdk \
    default-jre-headless \
    r-base \
    r-base-dev \
    gradle \
    cmake \
    make \
    gcc \
    g++ \
 && rm -rf /var/lib/apt/lists/*

export CFLAGS="-I/usr/include/python3.11 -I/usr/include/python3.11/cpython"
export CPPFLAGS="$CFLAGS"

# 2. Install Python packages
pip3 install --no-cache-dir --break-system-packages \
    cutadapt \
    pysam \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    deeptools==3.5.0 \
    ont-fast5-api==3.3.0 \
    h5py

# 3. Install bioinformatics binaries
## samtools & minimap2
# Prefer system versions if available; otherwise build from source
command -v samtools >/dev/null 2>&1 || { \
    git clone --depth 1 https://github.com/samtools/samtools.git /tmp/samtools \
    && cd /tmp/samtools \
    && autoheader && autoconf -Wno-syntax -Wno-error \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) && make install \
    && cd / && rm -rf /tmp/samtools; }
command -v minimap2 >/dev/null 2>&1 || { \
    git clone --depth 1 https://github.com/lh3/minimap2.git /tmp/minimap2 \
    && cd /tmp/minimap2 \
    && make -j$(nproc) \
    && cp minimap2 /usr/local/bin/ \
    && cp *.py /usr/local/bin/ \
    && cd / && rm -rf /tmp/minimap2; }

## seqkit
command -v seqkit >/dev/null 2>&1 || { \
    curl -Lo /usr/local/bin/seqkit \
         https://github.com/shenwei356/seqkit/releases/download/v0.11.0/seqkit_linux_amd64 \
    && chmod +x /usr/local/bin/seqkit; }

# 4. Install Jvarkit
if [ ! -f /usr/local/bin/jvarkit.jar ]; then
    git clone https://github.com/lindenb/jvarkit.git /tmp/jvarkit \
    && cd /tmp/jvarkit \
    && git checkout 014d3e9 \
    && ./gradlew biostar84452 \
    && cp dist/biostar84452.jar /usr/local/bin/jvarkit.jar \
    && cd / && rm -rf /tmp/jvarkit;
fi

# 5. Install Nanopolish
if [ ! -f /usr/local/bin/nanopolish ]; then
    git clone --recursive https://github.com/jts/nanopolish.git /tmp/nanopolish \
    && cd /tmp/nanopolish \
    && git checkout 480fc85 \
    && make -j$(nproc) \
    && cp nanopolish /usr/local/bin/ \
    && cd / && rm -rf /tmp/nanopolish;
fi

# 6. Install Perl modules (system-wide)
# Using cpanminus (cpanm)
cpanm --notest \
    Modern::Perl \
    Getopt::Long::Descriptive \
    Params::Validate \
    Params::Util \
    Sub::Install \
    autodie \
    Devel::Size \
    IO::File \
    IO::Interactive \
    IO::Uncompress::Gunzip \
    DBI \
    MooseX::App::Simple \
    MooseX::App::Command \
    MooseX::Getopt::Meta::Attribute::Trait::NoGetopt \
    CLIPSeqTools

# 7. Install R packages
Rscript -e 'install.packages(c("longitudinal", "fdrtool"), repos="http://cran.us.r-project.org")'
Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/GeneCycle/GeneCycle_1.1.5.tar.gz", repos=NULL, type="source")'

# 8. Install sam_to_sqlite and annotate-sqlite-with-fastq
# Assuming these scripts are in tools/utils
cd "${benchmark_dir}"
chmod +x tools/utils/sam_to_sqlite tools/utils/annotate-sqlite-with-fastq

# Cleanup
apt-get clean && rm -rf /var/lib/apt/lists/ /tmp/*
