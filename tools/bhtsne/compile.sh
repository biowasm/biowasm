#!/bin/bash

# Unzip sample data; we'll preload it with the module for convenience.
# -k keeps the .gz so the tracked file is never deleted/recreated (a re-gzip would rewrite
# the gzip mtime header and show a spurious diff even though the data is identical).
gunzip -kf ../data/brain8.snd.gz
gunzip -kf ../data/pollen2014.snd.gz

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
    CFLAGS="-O3 -s USE_ZLIB=1 -w" \
    LIBS="-s USE_ZLIB=1 -lm $FLAGS"

# Remove the extracted copies so the working tree stays clean (the .gz files were kept above).
rm -f ../data/brain8.snd ../data/pollen2014.snd
