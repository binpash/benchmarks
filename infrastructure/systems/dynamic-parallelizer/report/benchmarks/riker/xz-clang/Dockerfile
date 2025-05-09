FROM hs-riker

RUN apt-get update && apt-get install -y git \
    gcc \
    clang

WORKDIR /

RUN git clone https://github.com/xz-mirror/xz
WORKDIR /xz
RUN git checkout 2327a461e1afce862c22269b80d3517801103c1b

COPY xz-clang_build.sh /xz/xz-clang_build.sh
COPY wrapper.sh /xz/wrapper.sh
COPY xz-clang_build.sh /xz/Rikerfile
RUN chmod +x /xz/Rikerfile
RUN chmod +x /xz/wrapper.sh
RUN chmod +x /xz/xz-clang_build.sh