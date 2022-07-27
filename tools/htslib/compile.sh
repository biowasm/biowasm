#!/bin/bash

# TODO: look into LZMA support

# Dependencies
sudo apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev autoconf

# Run ./configure
make clean
autoheader
autoconf
emconfigure ./configure CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1" --disable-lzma

# Build
TOOLS=("tabix" "htsfile" "bgzip")
for tool in ${TOOLS[@]}; do
    emmake make $tool CC=emcc AR=emar \
        CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1" \
        LDFLAGS="$EM_FLAGS -O2"
done
