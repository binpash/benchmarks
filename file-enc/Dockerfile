FROM debian:12.7

WORKDIR /covid

COPY . .

RUN ./deps.sh
RUN ./input.sh
RUN ./run.sh
RUN ./verify.sh
RUN ./cleanup.sh
