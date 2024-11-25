IN=$PASH_TOP/aurpkg/input
mkdir -p ${IN}/deps/
# install dependencies
pkgs='ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump'

if ! dpkg -s $pkgs >/dev/null 2>&1 ; then
    sudo apt-get install $pkgs -y
    echo 'Packages Installed'
fi

# Check if makedeb-makepkg is installed
if ! dpkg -s "makedeb-makepkg" >/dev/null 2>&1 ; then
    cd ${IN}/deps/
    wget https://shlink.makedeb.org/install -O install.sh
    chmod +x install.sh
    ./install.sh
    echo 'Makedeb installed'
fi
