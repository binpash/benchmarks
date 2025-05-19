#!/usr/bin/env bash
set -euo pipefail

TOP=$(git rev-parse --show-toplevel)
benchmark_dir="${TOP}/teraseq"

# install.sh: Installs system-wide dependencies for the TERA-Seq pipeline

# 1. Install OS packages
apt-get update && sudo apt-get install -y --no-install-recommends \
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
    gffread \
    gmap \
 && rm -rf /var/lib/apt/lists/*

sudo wget -qO /usr/local/bin/liftOver http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
sudo chmod +x /usr/local/bin/liftOver

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
command -v seqkit >/dev/null 2>&1 || {
  curl -L https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz \
    -o /tmp/seqkit.tar.gz \
  && \
  # extract just the `seqkit` executable into /tmp
  tar -xzf /tmp/seqkit.tar.gz -C /tmp seqkit \
  && \
  mv /tmp/seqkit /usr/local/bin/seqkit \
  && chmod +x /usr/local/bin/seqkit \
  && rm /tmp/seqkit.tar.gz
}

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
    inc::Module::Install@1.19 \
    autodie@2.29 \
    Modern::Perl@1.20190601 \
    Getopt::Long::Descriptive@0.104 \
    Params::Validate@1.29 \
    Params::Util@1.07 \
    Sub::Install@0.928 \
    Devel::Size@0.83 \
    IO::File@1.39 \
    IO::Interactive@1.022 \
    IO::Uncompress::Gunzip \
    DBI@1.642 \
    MooseX::App::Simple@1.41 \
    MooseX::App::Command \
    MooseX::Getopt::Meta::Attribute::Trait::NoGetopt@0.74

cpanm --notest --force \
    GenOO@1.5.2 \
    CLIPSeqTools@0.1.9

# Manually install GenOOx from TeRA-Seq
tmpdir=$(mktemp -d) || { echo "Failed to create tempdir"; exit 1; }
git clone --depth 1 https://github.com/mourelatos-lab/TERA-Seq_manuscript.git "$tmpdir"
perldir=$(perl -MConfig -e 'print $Config{installsitelib}')
cp -r "$tmpdir"/misc/GenOOx "$perldir"
rm -rf "$tmpdir"

# 7. Install R packages
Rscript -e 'install.packages(c(
    "DBI",        
    "RSQLite",    
    "dplyr",      
    "stringr",    
    "optparse",   
    "longitudinal", 
    "fdrtool"     
  ), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/GeneCycle/GeneCycle_1.1.5.tar.gz", repos=NULL, type="source")'

# 8. Install sam_to_sqlite and annotate-sqlite-with-fastq
# Assuming these scripts are in tools/utils
cd "${benchmark_dir}"
chmod +x tools/utils/sam_to_sqlite tools/utils/annotate-sqlite-with-fastq

# Cleanup
apt-get clean && rm -rf /var/lib/apt/lists/ /tmp/*
