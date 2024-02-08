#!/bin/bash

make clean

# LZMA library flags
LZMA_VERSION="5.2.5"
DIR_LZMA=../../htslib/src/xz-${LZMA_VERSION}/src/liblzma
CFLAGS_LZMA="-I${DIR_LZMA}/api -I${DIR_LZMA}/api/lzma"
LFDLAGS_LZMA="-L${DIR_LZMA}/.libs"

# Use autoheader/autoconf to generate config.h.in and configure
autoheader
autoconf -Wno-syntax
emconfigure ./configure --without-curses --with-htslib="../../htslib/src/" CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"

# Build
emmake make samtools CC=emcc AR=emar \
    CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1 $CFLAGS_LZMA" \
    LDFLAGS="$EM_FLAGS --preload-file examples/@/samtools/examples/ -s ERROR_ON_UNDEFINED_SYMBOLS=0 -O2 $LFDLAGS_LZMA"
