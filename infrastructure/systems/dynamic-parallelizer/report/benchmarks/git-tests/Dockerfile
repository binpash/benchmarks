FROM hs
RUN apt update && apt-get install -q -y gettext libsvn-perl apache2 cvs cvsps git gnupg subversion make libssl-dev libcurl4-openssl-dev libexpat-dev wget sudo default-jre tcl tk gettext zlib1g-dev perl-modules liberror-perl libauthen-sasl-perl libemail-valid-perl libio-pty-perl libio-socket-ssl-perl libnet-smtp-ssl-perl libdbd-sqlite3-perl libcgi-pm-perl
RUN useradd -m hsuser
USER hsuser
WORKDIR /home/hsuser/
RUN git clone https://github.com/git/git.git
WORKDIR /home/hsuser/git
RUN git checkout bea9ecd24b0c3bf06cab4a851694fe09e7e51408 # v2.45.2
COPY --chown=hsuser:hsuser ./git_diff ./
RUN patch -p1 < ./git_diff
RUN make configure && ./configure && make -j $(nproc)
WORKDIR /home/hsuser/git/t
COPY --chown=hsuser:hsuser script_names.txt script_names.txt
COPY --chown=hsuser:hsuser run_tests.sh run_tests.sh

