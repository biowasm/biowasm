#!/bin/bash

# TODO: look into LZMA support

test -d ../htslib/build/ || (echo "Running 'make htslib' first." && cd ../../ && make htslib && cd tools/samtools/)

cd src/

make clean

# Also, use autoheader/autoconf to generate config.h.in and configure
autoheader
autoconf -Wno-syntax
emconfigure ./configure --without-curses --with-htslib="../../htslib/src/" CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"
emmake make CC=emcc AR=emar CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1"

# Rename output to .o so it's recognizable by Emscripten
cp samtools samtools.o

# Generate .wasm/.js files
emcc -O2 samtools.o \
    -o ../build/samtools.html \
    $EM_FLAGS \
    -s USE_BZIP2=1 \
    -s ERROR_ON_UNDEFINED_SYMBOLS=0 \
    --preload-file examples/@/tmp/examples/
