# Download base image ubuntu 16.04
FROM ubuntu:16.04

# LABEL about the custom image
LABEL maintainer="jan.oppelt@pennmedicine.upenn.edu"
LABEL version="0.3"
LABEL description="This is custom Docker Image for \
analysis of TERA-Seq publication (DOI: https://doi.org/10.1093/nar/gkab713)."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# Set default shell
SHELL ["/bin/bash", "-c"]

### System-wide requirements; cpanminus is not required if Perl uses virtual environment method; g++, zlib1g-dev, and bzip2 are required only for Nanopolish
RUN apt-get update \
    && apt-get install -y git gcc make wget g++ zlib1g-dev bzip2 \
    && rm -rf /var/lib/apt/lists/*

### Main GitHub repo
WORKDIR /root
COPY . TERA-Seq_manuscript

### Install Miniconda3
ENV PATH "/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O Miniconda3.sh \
    && mkdir /root/.conda \
    && bash Miniconda3.sh -b \
    && rm -f Miniconda3.sh
#RUN conda --version

## Install Mamba for faster installation
#RUN conda install -c conda-forge mamba

# Get Conda yml and install environment
#RUN mamba env create -f /usr/local/TERA-Seq_manuscript/teraseq-env.yml
RUN conda env create -f /root/TERA-Seq_manuscript/teraseq-env.yml

# Increase default FastQC RAM
RUN sed -i 's/-Xmx250m/-Xmx5g/g' /root/miniconda3/envs/teraseq/opt/fastqc-*/fastqc

#ENV PATH="${PATH}:/root/miniconda3/envs/teraseq/bin"

RUN ln -s /root/miniconda3/envs/teraseq/bin/R /bin/R \
    && ln -s /root/miniconda3/envs/teraseq/bin/curl /bin/curl

### Save default Conda path
RUN sed -i '/CONDA_PREFIX/d' /root/TERA-Seq_manuscript/PARAMS.sh \
    && echo -e "CONDA_PREFIX=\"/root/miniconda3\"" >> /root/TERA-Seq_manuscript/PARAMS.sh

### Perl
WORKDIR /root/TERA-Seq_manuscript/tools

## System-wide install (option 1)
# RUN cpanm inc::Module::Install \
#     && cpanm autodie \
#     && cpanm DBI \
#     && cpanm Devel::Size \
#     && cpanm Getopt::Long::Descriptive \
#     && cpanm IO::File \
#     && cpanm IO::Interactive \
#     && cpanm IO::Uncompress::Gunzip \
#     && cpanm Params::Validate \
#     && cpanm Params::Util \
#     && cpanm Sub::Install \
#     && cpanm Modern::Perl \
#     && cpanm --force MooseX::App::Simple \
#     && cpanm --force MooseX::App::Command \
#     && cpanm --force MooseX::Getopt::Meta::Attribute::Trait::NoGetopt
#
# WORKDIR /usr/local/share/perl/5.22.1
#
# RUN git clone --recursive https://github.com/genoo/GenOO.git GenOO_git
# WORKDIR /usr/local/share/perl/5.22.1/GenOO_git
# RUN git reset 6527029 --hard
# WORKDIR /usr/local/share/perl/5.22.1
# RUN mkdir GenOO \
#     && cp -r GenOO_git/lib/GenOO/* GenOO/
#
# RUN cpanm CLIPSeqTools
#
# RUN wget https://raw.githubusercontent.com/mourelatos-lab/TERA-Seq_manuscript/main/misc/GenOOx/Data/File/SAMminimap2.pm -O /usr/local/share/perl/5.22.1/GenOOx/Data/File/SAMminimap2.pm
# RUN mkdir /usr/local/share/perl/5.22.1/GenOOx/Data/File/SAMminimap2 \
#     && wget https://raw.githubusercontent.com/mourelatos-lab/TERA-Seq_manuscript/main/misc/GenOOx/Data/File/SAMminimap2/Record.pm -O /usr/local/share/perl/5.22.1/GenOOx/Data/File/SAMminimap2/Record.pm

## Virtual environment install (option 2)
# Export Conda perl lib path (mainly for local::lib module)
ENV PERL5LIB "/root/miniconda3/envs/teraseq/lib/site_perl/5.26.2/:${PERL5LIB}"
ARG PERL5LIB="/root/miniconda3/envs/teraseq/lib/site_perl/5.26.2/:${PERL5LIB}"

RUN git clone https://github.com/jizhang/perl-virtualenv.git \
    && cd perl-virtualenv/ \
    && git reset f931774 --hard \
    && chmod u+x virtualenv.pl \
    && ./virtualenv.pl teraseq \
    && . teraseq/bin/activate \
    && curl -L https://cpanmin.us/ -o teraseq/bin/cpanm \
    && chmod +x teraseq/bin/cpanm

RUN . perl-virtualenv/teraseq/bin/activate \
    && cpanm inc::Module::Install@1.19 \
    && cpanm autodie@2.29 \
    && cpanm DBI@1.642 \
    && cpanm Devel::Size@0.83 \
    && cpanm Getopt::Long::Descriptive@0.104 \
    && cpanm IO::File@1.39 \
    && cpanm IO::Interactive@1.022 \
    && cpanm --force IO::Uncompress::Gunzip \
    && cpanm Params::Validate@1.29 \
    && cpanm Params::Util@1.07 \
    && cpanm Sub::Install@0.928 \
    && cpanm Modern::Perl@1.20190601 \
    && cpanm --force MooseX::App::Simple@1.41 \
    && cpanm --force MooseX::App::Command \
    && cpanm --force MooseX::Getopt::Meta::Attribute::Trait::NoGetopt@0.74

RUN git clone --recursive https://github.com/genoo/GenOO.git perl-virtualenv/teraseq/lib/perl5/GenOO_git \
    && cd perl-virtualenv/teraseq/lib/perl5/GenOO_git/ \
    && git reset 6527029 --hard \
    && cd ../ \
    && mkdir GenOO \
    && cp -r GenOO_git/lib/GenOO/* GenOO/

# Install specific version of Perl module https://stackoverflow.com/questions/260593/how-can-i-install-a-specific-version-of-a-set-of-perl-modules
RUN . perl-virtualenv/teraseq/bin/activate \
    && cpanm --force CLIPSeqTools@0.1.9  \
    && cp -r /root/TERA-Seq_manuscript/misc/GenOOx/* perl-virtualenv/teraseq/lib/perl5/GenOOx/

################################################################################
### Nanopolish
# Default version
RUN git clone --recursive https://github.com/jts/nanopolish.git \
    && mv nanopolish nanopolish-480fc85 \
    && cd nanopolish-480fc85/ \
    && git reset 480fc85 --hard \
    && sed -i 's#http://bitbucket.org/eigen/eigen/get/$(EIGEN_VERSION).tar.bz2#https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.bz2#' Makefile \
    && sed -i 's/tar -xjf $(EIGEN_VERSION).tar.bz2/tar -xjf eigen-$(EIGEN_VERSION).tar.bz2/' Makefile \
    && sed -i 's/eigen-eigen-\*/eigen-$(EIGEN_VERSION)/' Makefile \
#    && sed -i '27 i EIGEN_VERSION_MV ?= d9c80169e091a2c6e75ceb509f81764d22cf6a63' Makefile \
#    && sed -i 's/mv\ eigen-\$(EIGEN_VERSION)/mv\ eigen-\$(EIGEN_VERSION_MV)/' Makefile \
    && rm -rf fast5 \
    && git clone https://github.com/mateidavid/fast5.git \
    && cd fast5/ \
    && git reset 18d6e34 --hard \
    && cd ../ \
    && rm -rf htslib \
    && git clone --recursive https://github.com/samtools/htslib.git \
    && cd htslib/ \
    && git reset 3dc96c5 --hard \
    && cd ../ \
    && make \
    && ln -s $(pwd)/nanopolish /root/miniconda3/envs/teraseq/bin/nanopolish

# New version with polya hmm scripts
RUN git clone --recursive https://github.com/jts/nanopolish.git \
    && mv nanopolish nanopolish-ab9722b \
    && cd nanopolish-ab9722b/ \
    && git reset ab9722b --hard

################################################################################
### Other dependencies
# Make sure to activate Conda
SHELL ["conda", "run", "-n", "teraseq", "/bin/bash", "-c"]

## GeneCycle
#RUN Rscript -e 'install.packages("GeneCycle", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages(c("longitudinal", "fdrtool"), repos = "http://cran.us.r-project.org"); install.packages("https://cran.r-project.org/src/contrib/GeneCycle_1.1.5.tar.gz", repos=NULL, type="source")'

## Cutadapt
RUN mkdir cutadapt-2.5 \
    && cd cutadapt-2.5/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && python3 -m pip install --upgrade pip \
    && pip3 install cutadapt==2.5 pysam numpy pandas matplotlib seaborn \
    && which cutadapt

## DeepTools
RUN mkdir deepTools-3.5.0 \
    && cd deepTools-3.5.0/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && python3 -m pip install --upgrade pip \
    && pip3 install wheel \
    && pip3 install deeptools==3.5.0 \
    && deeptools --version

## ONT-Fast5-API
RUN mkdir ont-fast5-api \
    && cd ont-fast5-api/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && pip install ont-fast5-api==3.3.0 h5py seaborn

## Jvarkit
RUN git clone "https://github.com/lindenb/jvarkit.git" \
    && mv jvarkit jvarkit-014d3e9 \
    && cd jvarkit-014d3e9/ \
    && git reset 014d3e9 --hard \
    && ./gradlew biostar84452 \
    && mkdir $CONDA_PREFIX/share/jvarkit \
    && ln -s $(pwd)/dist/biostar84452.jar /root/miniconda3/envs/teraseq/share/jvarkit/remove-softlip.jar


# Add utils dir to PATH
ENV PATH "/usr/local/TERA-Seq_manuscript/tools/utils:${PATH}"

WORKDIR /root/TERA-Seq_manuscript
