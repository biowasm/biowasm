#!/bin/bash

cd src/

# Unzip data; we'll preload it with the module for convenience
test -f brain8.snd || gunzip brain8.snd.gz

# Compile
rm t-sne.o
emmake make \
    CC=emcc CXX=em++ \
    CFLAGS+="-s USE_ZLIB=1" \
    LIBS="-s USE_ZLIB=1 -lm"

# Generate .wasm/.js files
mv t-sne t-sne.o
emcc -O2 -o ../build/t-sne.html t-sne.o \
    $EM_FLAGS \
    --preload-file brain8.snd@/bhtsne/brain8.snd
