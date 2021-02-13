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
    BT_LDFLAGS="--preload-file test@/bedtools2/test $EM_FLAGS"
