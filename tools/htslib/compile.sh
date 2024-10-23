#!/bin/bash

# Dependencies
sudo apt-get install -y zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev autoconf

# Compile LZMA to WebAssembly
LZMA_VERSION="5.6.3"
curl -LO "https://tukaani.org/xz/xz-${LZMA_VERSION}.tar.gz"
tar -xvf xz-${LZMA_VERSION}.tar.gz
cd xz-${LZMA_VERSION}
emconfigure ./configure --disable-shared --disable-threads
emmake make -j4 CFLAGS="-Oz -fPIC -s USE_PTHREADS=0 -s EXPORT_ALL=1 -s ASSERTIONS=1"
cd -

# Run ./configure
CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1 ${CFLAGS_LZMA}"
LDFLAGS="$LDFLAGS_LZMA"
make clean
autoheader
autoconf
emconfigure ./configure CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS"

# Build htslib tools
TOOLS=("tabix" "htsfile" "bgzip")
for tool in ${TOOLS[@]}; do
    emmake make $tool CC=emcc AR=emar \
        CFLAGS="-O2 $CFLAGS" \
        LDFLAGS="$EM_FLAGS -O2 -s ERROR_ON_UNDEFINED_SYMBOLS=0 $LDFLAGS"
done
