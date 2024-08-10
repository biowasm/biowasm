#!/bin/bash

# Download source from release to avoid needing yaggo and autotools
TARBALL=https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
mkdir -p src_release
curl -L -o mummer.tar.gz $TARBALL
tar -xzvf mummer.tar.gz -C src_release --strip-components=1

# Configure
cd src_release/
emconfigure ./configure --disable-openmp --disable-shared --host=none-none-none

# Remove multithreading
sed -i '152,156s/.*//' ./src/umd/nucmer_main.cc
sed -i 's/th.join();/query_thread(aligner.get(), \&parser, \&output, \&args);/' ./src/umd/nucmer_main.cc

# Compile nucmer
emmake make nucmer.js \
    EXEEXT=".js" \
    CFLAGS="-O2 $EM_FLAGS" \
    CXXFLAGS="-std=c++0x -O2 $EM_FLAGS"

mv nucmer.* ../../build
