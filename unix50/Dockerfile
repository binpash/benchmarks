FROM debian:12.7

WORKDIR /unix50

COPY . .

RUN apt-get update && apt-get install -y \
    bash \
    curl \
    python3 \
    python3-pip

RUN ./deps.sh
RUN ./input.sh
RUN ./run.sh
RUN ./verify.sh
RUN ./cleanup.sh

