## Running instructions

`git pull ghcr.io/binpash/hs/lua:latest`

Or build it yourself by running `docker build -t hs/lua .` from this directory.

Running with HS
`docker run --privileged -it --rm ghcr.io/binpash/hs/lua:latest /bin/bash -c 'time /srv/hs/pash-spec.sh /root/lua-5.4.3/lua_build.sh &> /dev/null && sha256sum /root/lua-5.4.3/src/lua'`

```
real    0m6.121s
user    0m5.379s
sys     0m1.111s
aeefc061ed3bb3f6831fcabf5e09f4f9dc53c5dbb91f0858230baf47bff25701  /root/lua-5.4.3/src/lua
```

Running with shell
`docker run --privileged -it --rm ghcr.io/binpash/hs/lua:latest /bin/bash -c 'time sh /root/lua-5.4.3/lua_build.sh &> /dev/null && sha256sum /root/lua-5.4.3/src/lua'`

```
real    0m3.155s
user    0m2.934s
sys     0m0.221s
aeefc061ed3bb3f6831fcabf5e09f4f9dc53c5dbb91f0858230baf47bff25701  /root/lua-5.4.3/src/lua
```
