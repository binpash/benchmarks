#/bin/bash

cd ..
apt install python3

function lsof_deps {
    python3 -m venv lsof_venv
    source lsof_venv/bin/activate
    apt install -y gcc libtirpc-dev
    export CPATH="$(pwd)/include:$(pwd)/lib/dialects/linux:$(pwd)/lib:/usr/include/tirpc:$(pwd)/src"
    deactivate
}

function lua_deps {
    python3 -m venv lua_venv
    source lua_venv/bin/activate
    apt install -y gcc make wget
    deactivate
}

function memcached_deps {
    python3 -m venv memcached_venv
    source memcached_venv/bin/activate
    apt install -y gcc autotools-dev automake libevent-dev
    deactivate
}

function redis_deps {
    python3 -m venv redis_venv
    source redis_venv/bin/activate
    apt install -y gcc make pkg-config tcl
    deactivate
}

function sqlite_deps {
    python3 -m venv sqlite_venv
    source sqlite_venv/bin/activate
    echo "Europe/London" > /etc/timezone
    dpkg-reconfigure -f noninteractive tzdata
    apt install -y gcc tcl8.6 libtool libtool-bin libreadline-dev tcl8.6-dev
    deactivate
}

function vim_deps {
    python3 -m venv vim_venv
    source vim_venv/bin/activate
    apt install -y gcc make libncurses-dev libsm-dev libice-dev libxt-dev libx11-dev libxdmcp-dev libselinux-dev
    deactivate
}

function xz_deps {
    python3 -m venv xz_venv
    source xz_venv/bin/activate
    apt install -y gcc
    deactivate
}

lsof_deps
lua_deps
memcached_deps
redis_deps
sqlite_deps
vim_deps
xz_deps