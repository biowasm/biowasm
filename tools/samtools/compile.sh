#!/bin/bash

# TODO: look into LZMA support

test -d ../htslib/build/ || echo "Run 'make htslib' first."

cd src/

make clean

# Also, use autoheader/autoconf to generate config.h.in and configure
autoheader
autoconf -Wno-syntax
emconfigure ./configure --without-curses --with-htslib="../../htslib/src/" CFLAGS="-s USE_ZLIB=1"
emmake make CC=emcc AR=emar

# Rename output to .o so it's recognizable by Emscripten
cp samtools samtools.o

# Generate .wasm/.js files 
emcc samtools.o \
    -o ../build/samtools.html \
    $EM_FLAGS \
    -s ERROR_ON_UNDEFINED_SYMBOLS=0 \
    --preload-file examples/@/tmp/examples/
