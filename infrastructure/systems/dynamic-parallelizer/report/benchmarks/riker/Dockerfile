FROM hs

WORKDIR /

RUN apt-get install -y make clang llvm git gcc python3-cram file graphviz
RUN git clone --recursive https://github.com/curtsinger-lab/riker.git

WORKDIR /riker
RUN make
RUN sudo make install