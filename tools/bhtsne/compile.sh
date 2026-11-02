#!/bin/bash

gunzip -c ../data/brain8.snd.gz > ../data/brain8.snd
gunzip -c ../data/pollen2014.snd.gz > ../data/pollen2014.snd

# Compile
FLAGS=$(cat <<EOF
    $EM_FLAGS \
    -O3 \
    -s ASYNCIFY=1 \
    -s 'ASYNCIFY_IMPORTS=["send_names","send_results"]' \
    -s EXPORTED_RUNTIME_METHODS=["getValue","UTF8ToString","callMain","FS","PROXYFS","WORKERFS","HEAPF64"] \
    --preload-file ../data/brain8.snd@/bhtsne/brain8.snd \
    --preload-file ../data/pollen2014.snd@/bhtsne/pollen2014.snd
EOF
)

emmake make \
    CC=emcc CXX=em++ \
    CFLAGS="-O2 -s USE_ZLIB=1 -w" \
    LIBS="-s USE_ZLIB=1 -lm $FLAGS"

# Drop the decompressed copies
cd -
rm -f data/brain8.snd data/pollen2014.snd
