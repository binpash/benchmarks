FROM debian:12.7

WORKDIR /max-temp

COPY . .

RUN ./deps.sh
RUN ./input.sh
RUN ./run.sh
RUN ./verify.sh
RUN ./cleanup.sh
