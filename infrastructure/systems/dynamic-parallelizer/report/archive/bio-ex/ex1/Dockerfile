FROM hs

WORKDIR /root
RUN mkdir Exercise1
WORKDIR /root/Exercise1

RUN apt-get update && apt-get install -y bedtools

RUN wget "http://rhodes.mm.di.uoa.gr:9999/DL/ex1/chr20.bed"
RUN wget "http://rhodes.mm.di.uoa.gr:9999/DL/ex1/chr20_flank500.bed"
RUN wget "http://rhodes.mm.di.uoa.gr:9999/DL/ex1/ATF3_chr20.bed"
RUN wget "http://rhodes.mm.di.uoa.gr:9999/DL/ex1/chr19.bed"
RUN wget "http://rhodes.mm.di.uoa.gr:9999/DL/ex1/ATF3_chr19.bed"

COPY ex1.sh /root/Exercise1/ex1.sh
RUN chmod +x /root/Exercise1/ex1.sh