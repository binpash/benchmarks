## Running instructions

`git pull ghcr.io/binpash/hs/redis:latest`

Or build it yourself by running `docker build -t hs/redis .` from this directory.

Running with HS
`docker run --privileged --rm ghcr.io/binpash/hs/redis:latest /bin/bash -c 'time /srv/hs/pash-spec.sh /root/redis/redis_build.sh && file /root/redis/src/redis-server'`

```
real    0m24.309s
user    0m20.368s
sys     0m3.903s
/root/redis/src/redis-server: ELF 64-bit LSB pie executable, x86-64, version 1 (SYSV), dynamically linked, interpreter /lib64/ld-linux-x86-64.so.2, BuildID[sha1]=36a1f18dbce7bf35eafe05da751f2d6906578341, for GNU/Linux 3.2.0, with debug_info, not stripped
```

Running with shell
`docker run --privileged --rm ghcr.io/binpash/hs/redis:latest /bin/bash -c 'time sh /root/redis/redis_build.sh && file /root/redis/src/redis-server'`

```
real    0m14.515s
user    0m13.439s
sys     0m0.800s
/root/redis/src/redis-server: ELF 64-bit LSB pie executable, x86-64, version 1 (SYSV), dynamically linked, interpreter /lib64/ld-linux-x86-64.so.2, BuildID[sha1]=e4bf8517ffe0c4edf32515ca936b6be5dc22a034, for GNU/Linux 3.2.0, with debug_info, not stripped
```
