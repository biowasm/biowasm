#!/bin/bash

cd src/

# Unzip data; we'll preload it with the module for convenience
test -f brain8.snd || gunzip brain8.snd.gz

# Compile to .wasm
emmake make \
    PROG="t-sne.html" \
    CC=emcc CXX=em++ \
    CFLAGS+="-s USE_ZLIB=1" \
    LIBS="-s USE_ZLIB=1 -lm --preload-file brain8.snd"

# Move files to build folder
for ext in data html js wasm; do
    mv "t-sne.${ext}" ../build/
done
