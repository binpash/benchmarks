# Running git test scripts

Building image with the Dockerfile `docker build -t hs/git-tests .`

Running image `docker run -it --rm --privileged --cgroupns=host hs/git-tests /bin/bash`
maybe also mount /tmp

installing gettext and other packages from https://github.com/git/git/blob/master/ci/install-dependencies.sh

t0211 has two failing cases, but I can't figure out how to get them successful.

Switched branch to release tag
