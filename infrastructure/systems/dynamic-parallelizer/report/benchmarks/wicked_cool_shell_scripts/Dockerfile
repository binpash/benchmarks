FROM hs

RUN apt update
RUN apt install -y time

WORKDIR /srv/hs
RUN mkdir -p report/benchmarks
COPY ./inner report/benchmarks/wicked_cool_shell_scripts
