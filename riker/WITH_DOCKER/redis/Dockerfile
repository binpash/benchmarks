FROM hs-riker

RUN apt update
RUN apt install -y git gcc make pkg-config tcl

WORKDIR /root
RUN git clone https://github.com/redis/redis

WORKDIR /root/redis
RUN git checkout d96f47cf06b1cc24b82109e0e87ac5428517525a
RUN make .make-prerequisites

COPY redis_build.sh /root/redis/redis_build.sh
COPY wrapper.sh /root/redis/wrapper.sh
COPY redis_build.sh /root/redis/Rikerfile
RUN chmod +x /root/redis/Rikerfile
RUN chmod +x /root/redis/wrapper.sh
RUN chmod +x /root/redis/redis_build.sh