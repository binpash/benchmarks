FROM hs

WORKDIR /root

RUN pip install --break-system-packages htseq
RUN apt-get install -y unzip default-jdk
RUN wget "https://uoa2-my.sharepoint.com/:u:/g/personal/nikosp_o365_uoa_gr/EUuyAqNvot9JlzRFJhmkq_8BQPbPM5z1-sb0CubekuDJEQ?e=tQg7ay&download=1" -O ex3.tar --no-check-certificate || wget http://rhodes.mm.di.uoa.gr:8888/ex3/ex3.tar -O ex3.tar
RUN tar -xvf ex3.tar

RUN git clone https://github.com/alexdobin/STAR.git
WORKDIR /root/STAR/source
RUN make

WORKDIR /root/workspace
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
RUN unzip Trimmomatic-0.39.zip
RUN unzip fastqc_v0.12.1.zip
RUN chmod +x /root/workspace/FastQC/fastqc

COPY modified_ex2.sh /root/workspace/modified_ex2.sh
RUN chmod +x /root/workspace/modified_ex2.sh