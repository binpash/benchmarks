## Running instructions

Build the image `docker build -t hs/calc .` from this directory.

Running with HS
`docker run --privileged --rm hs/calc:latest /bin/bash -c 'time /srv/hs/pash-spec.sh /root/calc/calc_build.sh && file /root/calc/calc'`

```
real    2m15.044s
user    2m9.292s
sys     0m52.000s
/root/calc/calc: ELF 64-bit LSB pie executable, x86-64, version 1 (SYSV), dynamically linked, interpreter /lib64/ld-linux-x86-64.so.2, BuildID[sha1]=e5bf6a86a5493060a783b6169950bba1265b4c70, for GNU/Linux 3.2.0, with debug_info, not stripped
```

Running with shell
`docker run --privileged --rm hs/calc:latest /bin/bash -c 'time sh /root/calc/calc_build.sh && file /root/calc/calc'`

```
real    0m8.590s
user    0m7.760s
sys     0m0.820s
/root/calc/calc: ELF 64-bit LSB pie executable, x86-64, version 1 (SYSV), dynamically linked, interpreter /lib64/ld-linux-x86-64.so.2, BuildID[sha1]=e5bf6a86a5493060a783b6169950bba1265b4c70, for GNU/Linux 3.2.0, with debug_info, not stripped
```
