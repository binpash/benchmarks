## Running instructions

Build the image `docker build -t hs/ex1 .` from this directory.

Running with HS
`docker run --privileged --rm hs/ex1:latest /bin/bash -c 'time /srv/hs/pash-spec.sh /root/Exercise1/ex1.sh && sha1sum /root/Exercise1/common_chr19.bed /root/Exercise1/common_chr20.bed'`

```
real    0m2.414s
user    0m2.705s
sys     0m2.340s
eb6b8df0b96a8eb0fe4cbad9bb8eef82e1dba09e  /root/Exercise1/common_chr19.bed
73e4af258f25ffd89d3e7e1621811b1cdc82dace  /root/Exercise1/common_chr20.bed
```

Running with shell
`docker run --privileged --rm hs/ex1:latest /bin/bash -c 'time sh /root/Exercise1/ex1.sh && sha1sum /root/Exercise1/common_chr19.bed /root/Exercise1/common_chr20.bed'`

```
real    0m0.037s
user    0m0.034s
sys     0m0.004s
eb6b8df0b96a8eb0fe4cbad9bb8eef82e1dba09e  /root/Exercise1/common_chr19.bed
73e4af258f25ffd89d3e7e1621811b1cdc82dace  /root/Exercise1/common_chr20.bed
```
