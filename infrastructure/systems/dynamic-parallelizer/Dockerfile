FROM debian:12

RUN mkdir -p /srv/hs
WORKDIR /srv/hs
SHELL ["/bin/bash", "-c"]
RUN apt update
RUN apt install -y vim sudo git python3 python3.11-venv strace wget make python3-cram file graphviz libtool python3-matplotlib libcap2-bin util-linux
# pash distro deps
RUN apt install -y bc curl graphviz bsdmainutils libffi-dev locales locales-all netcat-openbsd pkg-config procps python3-pip python3-setuptools python3-testresources wamerican-insane
# try deps
RUN apt install -y expect mergerfs attr
RUN git config --global --add safe.directory /srv
COPY . .
RUN python3 -m venv .venv
RUN source .venv/bin/activate
ENV PASH_SPEC_TOP=/srv/hs
ENV PASH_TOP=/srv/hs/deps/pash
WORKDIR /srv/hs/deps/try
RUN make -C utils
RUN mv utils/try-commit /bin
RUN mv utils/try-summary /bin
WORKDIR /srv/hs/deps/pash
RUN ./scripts/setup-pash.sh
WORKDIR /srv/hs
RUN chmod +x entrypoint.sh
ENTRYPOINT ["/srv/hs/entrypoint.sh"]
