#!/bin/sh

./autogen.sh
if [ ! -f Makefile ]; then
  ./configure
fi

make clean all