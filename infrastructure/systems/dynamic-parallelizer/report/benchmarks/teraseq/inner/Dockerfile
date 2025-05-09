FROM ghcr.io/binpash/teraseq20:latest
WORKDIR /root/TERA-Seq_manuscript
COPY . .
RUN mv -f activate.d/* /root/miniconda3/envs/teraseq/etc/conda/activate.d
RUN mv -f deactivate.d/* /root/miniconda3/envs/teraseq/etc/conda/deactivate.d
RUN ./inputs.sh
WORKDIR /root/TERA-Seq_manuscript/data
RUN ./run.sh
WORKDIR /root/TERA-Seq_manuscript/samples
