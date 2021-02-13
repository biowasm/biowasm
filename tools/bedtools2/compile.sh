#!/bin/bash

cd src/

# Remove large files (we'll pre-load the rest of the files as examples)
rm -rf ./test/intersect/sortAndNaming/bigTests

# Install dependencies
sudo apt-get install -y libbz2-dev liblzma-dev

make clean
# Compile to WebAssembly
emmake make \
    BIN_DIR="../build/" \
    CXX="em++ -s USE_ZLIB=1 -s USE_BZIP2=1" \
    BT_LIBS="-s USE_ZLIB=1 -s USE_BZIP2=1 -lm -lpthread" \
    BT_LDFLAGS="--preload-file test@/bedtools2/test $EM_FLAGS"
