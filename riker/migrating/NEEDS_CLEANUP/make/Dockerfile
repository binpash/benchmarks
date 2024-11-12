FROM hs-riker

RUN apt-get update && apt-get install -y wget \
    gcc \
    make

WORKDIR /

RUN wget http://ftp.gnu.org/gnu/make/make-4.3.tar.gz && \
    tar xzf make-4.3.tar.gz && \
    rm make-4.3.tar.gz
WORKDIR /make-4.3
RUN ./configure && \
    make clean

COPY make_build.sh /make-4.3/make_build.sh
COPY make_build.sh /make-4.3/Rikerfile
RUN chmod +x /make-4.3/Rikerfile
RUN chmod +x /make-4.3/make_build.sh