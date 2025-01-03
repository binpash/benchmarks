REPO_TOP=$(git rev-parse --show-toplevel)
IN=$REPO_TOP/aurpkg/input

mkdir -p ${IN}/deps/
pkgs='ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump makedeb'

# Add the makedeb repository
apt-get install gpg
wget -qO - 'https://proget.makedeb.org/debian-feeds/makedeb.pub' | gpg --dearmor | sudo tee /usr/share/keyrings/makedeb-archive-keyring.gpg 1> /dev/null
echo 'deb [signed-by=/usr/share/keyrings/makedeb-archive-keyring.gpg arch=all] https://proget.makedeb.org/ makedeb main' | sudo tee /etc/apt/sources.list.d/makedeb.list
sudo apt update

if ! dpkg -s $pkgs >/dev/null 2>&1 ; then
    sudo apt-get install $pkgs -y
    echo 'Packages Installed'
fi

useradd -m user && \
  echo "user ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
  chown -R user:user /benchmarks
