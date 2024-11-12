FROM hs-riker

RUN apt-get update && apt-get install -y git wget gcc make libgmp-dev perl libacl1-dev libattr1-dev libcap-dev

WORKDIR /

RUN wget --no-clobber https://ftp.gnu.org/gnu/coreutils/coreutils-8.32.tar.gz
RUN tar xzf coreutils-8.32.tar.gz
RUN rm coreutils-8.32.tar.gz

WORKDIR /coreutils-8.32
ENV FORCE_UNSAFE_CONFIGURE=1

RUN test -f Makefile || ./configure
COPY coreutils_build.sh /coreutils-8.32/coreutils_build.sh
COPY coreutils_build.sh /coreutils-8.32/Rikerfile
RUN chmod +x /coreutils-8.32/Rikerfile
RUN chmod +x /coreutils-8.32/coreutils_build.sh