#!/bin/bash

# Unzip sample data; we'll preload it with the module for convenience
test -f data/brain8.snd || gunzip data/brain8.snd.gz
test -f data/pollen2014.snd || gunzip data/pollen2014.snd.gz

# Compile
cd src/
rm t-sne.o
emmake make \
    CC=emcc CXX=em++ \
    CFLAGS+="-s USE_ZLIB=1" \
    LIBS="-s USE_ZLIB=1 -lm"

# Generate .wasm/.js files
mv t-sne t-sne.o
emcc -O2 -o ../build/t-sne.html t-sne.o \
    $EM_FLAGS \
    -s ASYNCIFY=1 \
    -s 'ASYNCIFY_IMPORTS=["send_names","send_results"]' \
    -s EXTRA_EXPORTED_RUNTIME_METHODS=["callMain","getValue","UTF8ToString"] \
    --preload-file ../data/brain8.snd@/bhtsne/brain8.snd \
    --preload-file ../data/pollen2014.snd@/bhtsne/pollen2014.snd

# Undo sample data unzipping
cd -
test -f data/brain8.snd.gz || gzip data/brain8.snd
test -f data/pollen2014.snd.gz || gzip data/pollen2014.snd
