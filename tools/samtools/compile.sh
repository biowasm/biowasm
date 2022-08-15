#!/bin/bash

# TODO: look into LZMA support

make clean

# Also, use autoheader/autoconf to generate config.h.in and configure
autoheader
autoconf -Wno-syntax
emconfigure ./configure --without-curses --with-htslib="../../htslib/src/" CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"

# Build
emmake make samtools CC=emcc AR=emar \
    CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1" \
    LDFLAGS="$EM_FLAGS --preload-file examples/@/samtools/examples/ -s ERROR_ON_UNDEFINED_SYMBOLS=0 -O2"
