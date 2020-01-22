#!/bin/bash

# TODO: look into LZMA support

# Dependencies
apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

# Run ./configure
cd src/
make clean
autoheader
autoconf
emconfigure ./configure CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1" --disable-lzma
