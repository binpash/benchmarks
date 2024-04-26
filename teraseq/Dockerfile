FROM ghcr.io/ezrizhu/teraseq:latest
VOLUME ["/root/TERA-Seq_manuscript/samples"]
WORKDIR /root/TERA-Seq_manuscript
COPY . .
WORKDIR /root/TERA-Seq_manuscript/data
RUN ./run.sh
