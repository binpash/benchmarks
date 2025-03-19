FROM debian:12.7

WORKDIR /log-analysis

COPY . .

RUN ./deps.sh
RUN ./input.sh
RUN ./run.sh
RUN ./verify.sh
RUN ./cleanup.sh
