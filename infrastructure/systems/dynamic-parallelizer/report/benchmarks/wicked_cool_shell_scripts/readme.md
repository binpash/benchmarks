* 27: alternative to cat -n
* 28: apply fmt to lines thats long
* 29: show the contents of a file, additional userful info
* 35: disk quota analysis for users - currently hangs
* 37: fancier version of df - output different due to df in try
* 40: adds a new user - fails, shadow file owned by another gid
* 41: suspends user - not tested, interactive focused
* 42: delete unix user - fails, shadow file owned by another gid
* 43: verify PATH and other env vars
* 48: verify cron file
* 50: log rotator
* 51: backup dirs - removal, timestamp involved, feature overlap
* 52: wrapper over tar czf
* 102: mass rename of files

`docker run --rm -it --privileged $(docker build -q .) /bin/bash -c "window=8 sh run.sh && bash"`
