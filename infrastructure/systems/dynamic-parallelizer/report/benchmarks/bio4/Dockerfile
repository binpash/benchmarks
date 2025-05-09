from hs

WORKDIR /srv/hs
RUN mkdir -p report/benchmarks/bio4
COPY ./inner/setup report/benchmarks/bio4/setup 
COPY ./inner/input_full.txt report/benchmarks/bio4/input_full.txt
RUN ./report/benchmarks/bio4/setup
COPY ./inner report/benchmarks/bio4
