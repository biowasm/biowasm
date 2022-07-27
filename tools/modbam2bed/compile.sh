#!/bin/bash

# Cribbed from htslib example
# we're using a development branch of htslib, so cannot use the biowasm version

sudo apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev autoconf

make clean

# build htslib.a
cd htslib
make clean
autoheader
autoconf
emconfigure ./configure CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1 -O2" --disable-lzma  --disable-libcurl
emmake make libhts.a CC="emcc" AR=emar
cd ..


# argp - this is a little gnarly
wget https://www.lysator.liu.se/~nisse/misc/argp-standalone-1.3.tar.gz
tar -xzvf argp-standalone-1.3.tar.gz
cd argp-standalone-1.3
wget https://raw.githubusercontent.com/Homebrew/formula-patches/b5f0ad3/argp-standalone/patch-argp-fmtstream.h
patch argp-fmtstream.h patch-argp-fmtstream.h
emconfigure ./configure 
emmake make CC=emcc AR=emar \
    LDFLAGS="$EM_FLAGS -O2"
cd ..


# build program
# The define `WASM` disables posix resource limit checks
# The define NOTHREADS disables use of hts_thread_pool which leads to error pthread_attr_getstacksize
# being missing when compilation is not performed without `-pthread`. The use of `-pthread` triggers
# warnings with emcc that pthreads + MODULARIZE requires setting EXPORT_NAME != Module. However this
# doesnt appear compatible with aioli. Consequently we need to set ERROR_ON_UNDEFINED_SYMBOLS=0 to
# ignore the undefined pthread_attr_getstacksize function.
emmake make modbam2bed CC=emcc AR=emar \
    CFLAGS="-DNOTHREADS -DWASM -s ERROR_ON_UNDEFINED_SYMBOLS=0 -s USE_ZLIB=1 -s USE_BZIP2=1 -fno-stack-protector -O2" \
    ARGP=argp-standalone-1.3/libargp.a \
    EXTRA_CFLAGS="-Iargp-standalone-1.3" \
    EXTRA_LDFLAGS="$EM_FLAGS"
cd ..
cp src/modbam2bed build/modbam2bed.js
cp src/modbam2bed.wasm build
