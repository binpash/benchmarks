apt update
apt install -y wget binutils gzip
mkdir /work
cd /work
mv /debootstrap_1.0.134_all.deb .
ar -x debootstrap_1.0.134_all.deb
cd /
zcat /work/data.tar.gz | tar xv
mkdir /srv/debinstall
