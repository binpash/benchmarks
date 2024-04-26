FROM hs-riker

RUN apt-get update
RUN apt-get install -y git gcc-11 make libreadline-dev bsdmainutils
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 1

WORKDIR /root
RUN git clone https://github.com/lcn2/calc.git

WORKDIR /root/calc
RUN git checkout 3d33c6c6f498d4d391918d34c55fe37191db53d8

COPY calc_build.sh /root/calc/calc_build.sh
COPY calc_build.sh /root/calc/Rikerfile
RUN chmod +x /root/calc/Rikerfile
RUN chmod +x /root/calc/calc_build.sh
