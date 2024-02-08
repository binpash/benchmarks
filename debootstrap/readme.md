## To build
`docker build -t debootstrap .`

## To run
Make the work dir for debootstrap to bootstrap debian.
`mkdir work`

Start the conainer to build debian.
`time docker run -it --rm -v ./work:/srv/debinstall debootstrap /bin/bash /download.sh`
`time docker run -it --rm -v ./work:/srv/debinstall debootstrap /bin/bash /install.sh`
