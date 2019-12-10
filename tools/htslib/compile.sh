#!/bin/bash

# TODO: figure out "--disable-bz2 --disable-lzma"

# Dependencies
apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

# Run ./configure
cd src/
autoheader
autoconf
emconfigure ./configure CFLAGS="-s USE_ZLIB=1" --disable-bz2 --disable-lzma
