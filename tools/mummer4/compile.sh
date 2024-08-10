#!/bin/bash

BRANCH_OR_TAG=$(git symbolic-ref -q --short HEAD || git describe --tags --exact-match)
TARBALL=https://github.com/mummer4/mummer/releases/download/v$BRANCH_OR_TAG/mummer-$BRANCH_OR_TAG.tar.gz

# Download source from release to avoid needing yaggo and autotools
mkdir -p src_release
wget -O mummer.tar.gz $TARBALL
tar -xzvf mummer.tar.gz -C src_release --strip-components=1

# Configure
cd src_release/
emconfigure ./configure --disable-openmp --disable-shared

# Remove multithreading
sed -i '152,156s/.*//' ./src/umd/nucmer_main.cc
sed -i 's/th.join();/query_thread(aligner.get(), \&parser, \&output, \&args);/' ./src/umd/nucmer_main.cc

# Compile nucmer
emmake make nucmer.js \
    EXEEXT=".js" \
    CFLAGS="-O2 $EM_FLAGS" \
    CXXFLAGS="-std=c++0x -O2 $EM_FLAGS"

mv nucmer.* ../../build
