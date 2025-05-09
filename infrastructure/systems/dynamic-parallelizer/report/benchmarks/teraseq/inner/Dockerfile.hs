FROM teraseq20-data
RUN mkdir -p /srv/hs
WORKDIR /srv/hs
SHELL ["/bin/bash", "-c"]
# ARG DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN apt update
RUN apt install -y sudo curl git
RUN conda install -y -c anaconda python=3.11
RUN apt install -y vim strace make python3-cram file graphviz libtool python3-matplotlib libcap2-bin mergerfs lsb-release sudo git python3 wget util-linux
# pash distro deps
RUN apt install -y bc curl graphviz bsdmainutils libffi-dev locales locales-all netcat-openbsd pkg-config procps python3-pip python3-setuptools python3-testresources wamerican-insane
# try deps
RUN apt install -y expect mergerfs attr
COPY deps deps
COPY parallel-orch parallel-orch
COPY .git .git
RUN mkdir -p /srv/hs/report/benchmarks/teraseq
COPY report/benchmarks/teraseq/inner /srv/hs/report/benchmarks/teraseq
RUN /srv/hs/report/benchmarks/teraseq/setup

RUN git config --global --add safe.directory /srv
RUN git submodule update --init --recursive
ENV PASH_SPEC_TOP=/srv/hs
ENV PASH_TOP=/srv/hs/deps/pash
WORKDIR /srv/hs/deps/try
RUN make -C utils
RUN mv utils/try-commit /bin
RUN mv utils/try-summary /bin
WORKDIR /srv/hs/deps/pash
RUN ./scripts/setup-pash.sh

WORKDIR /srv/hs
COPY scripts scripts
COPY entrypoint.sh entrypoint.sh
COPY pash-spec.sh pash-spec.sh
RUN chmod +x entrypoint.sh

RUN mv /srv/hs/report/benchmarks/teraseq/annotate-sqlite-with-fastq.R /root/TERA-Seq_manuscript/tools/utils/

