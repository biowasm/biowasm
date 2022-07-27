#!/bin/bash

# Unzip sample data; we'll preload it with the module for convenience
test -f ../data/brain8.snd || gunzip ../data/brain8.snd.gz
test -f ../data/pollen2014.snd || gunzip ../data/pollen2014.snd.gz

# Compile
FLAGS=$(cat <<EOF
    $EM_FLAGS \
    -O3 \
    -s ASYNCIFY=1 \
    -s 'ASYNCIFY_IMPORTS=["send_names","send_results"]' \
    -s EXPORTED_RUNTIME_METHODS=["getValue","UTF8ToString","callMain","FS","PROXYFS","WORKERFS"] \
    --preload-file ../data/brain8.snd@/bhtsne/brain8.snd \
    --preload-file ../data/pollen2014.snd@/bhtsne/pollen2014.snd
EOF
)

emmake make \
    CC=emcc CXX=em++ \
    CFLAGS="-O2 -s USE_ZLIB=1 -w" \
    LIBS="-s USE_ZLIB=1 -lm $FLAGS"

# Undo sample data unzipping
cd -
test -f data/brain8.snd.gz || gzip data/brain8.snd
test -f data/pollen2014.snd.gz || gzip data/pollen2014.snd
