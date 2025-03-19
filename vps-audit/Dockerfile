FROM debian:12.7

WORKDIR /vps-audit

COPY . .

RUN ./deps.sh
RUN ./input.sh
RUN ./run.sh
RUN ./verify.sh
RUN ./cleanup.sh
