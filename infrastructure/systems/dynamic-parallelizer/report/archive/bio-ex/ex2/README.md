## Running instructions

Build the image `docker build -t hs/ex2 .` from this directory.

Running with HS
`docker run --privileged --rm hs/ex2:latest /bin/bash -c 'time /srv/hs/pash-spec.sh /root/workspace/ex2.sh && sha1sum /root/workspace/htseq_output/SRR10045016-17-18-19-20-21_counts.csv'`

```
real    57m53.638s
user    120m52.257s
sys     5m53.287s
ded1b425f7b5141336b8759dc8a3b5eb68075269  /root/workspace/htseq_output/SRR10045016-17-18-19-20-21_counts.csv
```

Running with shell
`docker run --privileged --rm hs/ex2:latest /bin/bash -c 'time sh /root/workspace/ex2.sh && sha1sum /root/workspace/htseq_output/SRR10045016-17-18-19-20-21_counts.csv'`

```
real    48m49.509s
user    64m56.203s
sys     0m59.038s
ded1b425f7b5141336b8759dc8a3b5eb68075269  htseq_output/SRR10045016-17-18-19-20-21_counts.csv
```

hs with script modification

```
real    24m58.410s
user    129m29.244s
sys     5m58.171s
57771e0e4e2fd6b9d8c8a057bfe3948f89e83dee  /root/workspace/htseq_output/SRR10045016-17-18-19-20-21_counts.csv
```