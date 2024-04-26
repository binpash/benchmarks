#/bin/bash

cd ..
git clone https://github.com/lsof-org/lsof.git
cd lsof
git checkout 005e014e1abdadb2493d8b3ce87b37a2c0a2351d
./Configure -n linux

cd ..
wget https://www.lua.org/ftp/lua-5.4.3.tar.gz
tar xzf lua-5.4.3.tar.gz && mv lua-5.4.3 lua

git clone https://github.com/memcached/memcached.git
cd memcached
git checkout c472369fed5981ba8c004d426cee62d5165c47ca
./autogen.sh
./configure --disable-dependency-tracking

cd ..
git clone https://github.com/redis/redis
cd redis
git checkout d96f47cf06b1cc24b82109e0e87ac5428517525a
make .make-prerequisites

cd ..
git clone https://github.com/sqlite/sqlite
cd sqlite
git checkout c1cace0832fa2af5ab8315e217d708c09d586425
./configure

cd ..
git clone https://github.com/vim/vim.git
cd vim 
git checkout b836f631dba2534efd314a8f77439cebc75acd4e
./configure

cd ..
git clone https://github.com/xz-mirror/xz
cd xz
git checkout 2327a461e1afce862c22269b80d3517801103c1b