FROM hs-riker

RUN apt update
RUN apt install -y gcc make wget

WORKDIR /root

RUN wget https://www.lua.org/ftp/lua-5.4.3.tar.gz
RUN tar xzf lua-5.4.3.tar.gz

COPY lua_build.sh lua-5.4.3/lua_build.sh
COPY wrapper.sh lua-5.4.3/wrapper.sh
COPY lua_build.sh lua-5.4.3/Rikerfile
RUN chmod +x lua-5.4.3/Rikerfile
RUN chmod +x lua-5.4.3/lua_build.sh
RUN chmod +x lua-5.4.3/wrapper.sh

WORKDIR /root/lua-5.4.3
