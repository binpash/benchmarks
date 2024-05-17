#!/bin/bash

cd ..
apt install python3

lsof_deps () {
    echo "lsof"
    apt install -y gcc libtirpc-dev
    export CPATH_TEMP=$CPATH
    export CPATH="$(pwd)/include:$(pwd)/lib/dialects/linux:$(pwd)/lib:/usr/include/tirpc:$(pwd)/src"
}

lua_deps () {
    echo "lua"
    apt install -y gcc make wget
}

memcached_deps () {
    echo "memcached"
    apt install -y gcc autotools-dev automake libevent-dev
}

redis_deps () {
    echo "redis"
    apt install -y gcc make pkg-config tcl
}

sqlite_deps () {
    echo "sqlite"
    apt install -y gcc tcl8.6 libtool libtool-bin libreadline-dev tcl8.6-dev
}

vim_deps () {
    echo "vim"
    apt install -y gcc make libncurses-dev libsm-dev libice-dev libxt-dev libx11-dev libxdmcp-dev libselinux-dev
}

xz_deps () {
    echo "xz"
    apt install -y gcc
}

lsof_deps
lua_deps
memcached_deps
redis_deps
sqlite_deps
vim_deps
xz_deps