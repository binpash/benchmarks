FROM debian:12

COPY debootstrap_1.0.134_all.deb .
COPY prepare.sh .
COPY download.sh .
COPY install.sh .
RUN /bin/bash prepare.sh
